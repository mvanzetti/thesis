from pybedtools import BedTool
import random
from datetime import datetime
import pandas as pd

import overlaps_db.constants.dataframe_info as info


class OverlapAnalyzer:
    # TODO duplicate of overlapper.py method
    def compute_size(self, row, prefix=None):
        col_name = prefix + '_name' if prefix else 'name'
        if row[col_name] == '.':
            return 0
        col_end = prefix + '_end' if prefix else 'end'
        col_start = prefix + '_start' if prefix else 'start'
        size = abs(row[col_end] - row[col_start])
        return size

    def mean_size(self, bed):
        sample_df = bed.to_dataframe()
        sample_df['size'] = sample_df.apply(lambda row: self.compute_size(row), axis=1)
        return int(round(sample_df[['size']].mean()))

    def intersect(self, bed, bed_overlap_with, left_outer_join, min_overlap):
        return bed.intersect(bed_overlap_with, loj=left_outer_join, f=min_overlap)

    def create_null_overlap_model_by_shuffle(self, bed, bed_overlap_with, assembly, left_outer_join, min_overlap):
        print("Creating Null Model By Shuffle Strategy...")

        random_bed = bed_overlap_with.shuffle(genome=assembly, chrom=False, seed=random.seed(datetime.now()))

        random_overlap = self.intersect(bed, random_bed, left_outer_join, min_overlap)
        return random_overlap

    def create_null_overlap_model_by_random(self, bed, bed_overlap_with, assembly, left_outer_join, min_overlap):
        print("Creating Null Model By Random Strategy...")

        df_overlap_with = bed_overlap_with.to_dataframe()
        df_overlap_with['size'] = df_overlap_with.apply(lambda row: self.compute_size(row), axis=1)
        interval_size = int(round(df_overlap_with[['size']].mean()))
        interval_num = len(df_overlap_with)

        print("Mean interval size is", interval_size)
        print("Intervals are", interval_num)

        empty_bed = BedTool()
        random_bed = empty_bed.random(l=interval_size, n=interval_num, genome=assembly)

        random_overlap = self.intersect(bed, random_bed, left_outer_join, min_overlap)
        return random_overlap

    # create random sequence and overlap n times and compute a statistics
    # TODO parallelize
    def create_random_overlap_distribution(self, bed, bed_overlap_with, assembly, min_overlap, samples_num=20,
                                           strategy='random'):
        columns = ['sample_num', 'size']
        tests_df = pd.DataFrame(columns=columns)
        for i in range(0, samples_num):

            if strategy == 'shuffle':
                random_bed = bed_overlap_with.shuffle(genome=assembly, chrom=False)

            elif strategy == 'random':
                empty_bed = BedTool()
                random_bed = empty_bed.random(l=self.mean_size(bed_overlap_with), n=bed_overlap_with.count(),
                                              genome=assembly)
            else:
                raise ValueError('The specified strategy does not exist:', strategy)

            row_array = [i, bed.intersect(random_bed, f=min_overlap).count()]
            temp_df = pd.DataFrame([row_array], columns=columns)
            tests_df = tests_df.append(temp_df)
        return tests_df

    def compute_reldist(self, bed, bed_overlap_with):
        # TODO implement
        pass

    def compute_fisher_jaccard_tests(self, bed, bed_overlap_with, bed_name, bed_overlap_with_name, biosample_name,
                                     assembly, samples_num=10):
        columns = info.overlap_tests['fisher_jaccard'].columns

        tests_df = pd.DataFrame(columns=columns)

        a_size = len(bed)
        b_size = len(bed_overlap_with)

        for i in range(1, samples_num):
            min_ovlp = i * 1. / samples_num

            # fisher test
            fisher = bed.fisher(bed_overlap_with, f=min_ovlp, genome=assembly)
            overlaps_count = fisher.table['in -a']['in -b']
            left_tail_fisher_pvalue = fisher.left_tail
            right_tail_fisher_pvalue = fisher.right_tail
            two_tail_fisher_pvalue = fisher.two_tail
            oddsratio_fisher = fisher.ratio

            # jaccard index
            jaccard = bed.jaccard(bed_overlap_with, f=min_ovlp)
            jaccard_index = jaccard['jaccard']

            row_array = [bed_name, biosample_name, bed_overlap_with_name, a_size, b_size,
                         min_ovlp, overlaps_count, left_tail_fisher_pvalue, right_tail_fisher_pvalue,
                         two_tail_fisher_pvalue, oddsratio_fisher, jaccard_index]

            temp_df = pd.DataFrame([row_array], columns=columns)
            tests_df = tests_df.append(temp_df)

        tests_df.reset_index(inplace=True, drop=True)
        return tests_df

    def compute_fisher_jaccard_z_tests(self, bed, bed_overlap_with, bed_name, bed_overlap_with_name, biosample_name,
                                       assembly, samples_num=10, random_samples_num=20):
        columns = info.overlap_tests['fisher_jaccard_z'].columns

        tests_df = pd.DataFrame(columns=columns)

        a_size = len(bed)
        b_size = len(bed_overlap_with)

        for i in range(1, samples_num):
            min_ovlp = i * 1. / samples_num

            # fisher test
            fisher = bed.fisher(bed_overlap_with, f=min_ovlp, genome=assembly)
            overlaps_count = fisher.table['in -a']['in -b']
            right_tail_fisher_pvalue = fisher.right_tail

            # jaccard index
            jaccard = bed.jaccard(bed_overlap_with, f=min_ovlp)
            jaccard_index = jaccard['jaccard']

            # z test
            runs_shuffled_df = self.create_random_overlap_distribution(
                bed, bed_overlap_with, assembly, min_ovlp, random_samples_num, 'shuffle')
            runs_random_df = self.create_random_overlap_distribution(
                bed, bed_overlap_with, assembly, min_ovlp, random_samples_num)

            random_mean_count = float(runs_random_df[['size']].mean())
            random_std = float(runs_random_df[['size']].std())
            shuffled_mean_count = float(runs_shuffled_df[['size']].mean())
            shuffled_std = float(runs_shuffled_df[['size']].std())

            z_random = (overlaps_count - random_mean_count) / random_std
            z_shuffled = (overlaps_count - shuffled_mean_count) / shuffled_std

            row_array = [bed_name, biosample_name, bed_overlap_with_name, a_size, b_size,
                         min_ovlp, overlaps_count, z_random, z_shuffled, right_tail_fisher_pvalue, jaccard_index]

            temp_df = pd.DataFrame([row_array], columns=columns)
            tests_df = tests_df.append(temp_df)

        tests_df.reset_index(inplace=True, drop=True)
        return tests_df
