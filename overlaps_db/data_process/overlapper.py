from pybedtools import BedTool
import random
from datetime import datetime
import pandas as pd


class Overlapper:
    def __init__(self, download_path, staging_path, overlap_path):
        self.download_path = download_path
        self.staging_path = staging_path
        self.overlap_path = overlap_path

    @staticmethod
    def set_col_value(row, col_name, value):
        if row[col_name] == '.':
            return '.'
        return value

    @staticmethod
    def compute_size(row, prefix=None):
        col_name = prefix + '_name' if prefix else 'name'
        if row[col_name] == '.':
            return 0
        col_end = prefix + '_end' if prefix else 'end'
        col_start = prefix + '_start' if prefix else 'start'
        size = abs(row[col_end] - row[col_start])
        return size

    @staticmethod
    def compute_ovlp_len(row, prefix):
        col_name = prefix + '_name'
        col_end = prefix + '_end'
        col_start = prefix + '_start'
        if row[col_name] == '.':
            return 0
        value = min(row['end'], row[col_end]) - max(row['start'], row[col_start])
        return value

    @staticmethod
    def compute_ovlp_pct(row, prefix):
        col_len = prefix + '_ovlp_len'
        size = abs(row['end'] - row['start'])
        value = row[col_len] / size * 100.0
        return value

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

    def build_tests_df(self, bed, bed_overlap_with, bed_encyclopedia, bed_overlap_with_encyclopedia, biosample_name,
                       assembly):
        columns = ['encyclopedia', 'biosample_name', 'ovlp_encyclopedia', 'encyclopedia_size', 'ovlp_encyclopedia_size',
                   'min_ovlp', 'ovlp_count', 'fisher_left_p', 'fisher_right_p', 'fisher_two_p', 'jaccard']

        tests_df = pd.DataFrame(columns=columns)

        a_size = len(bed)
        b_size = len(bed_overlap_with)

        for i in range(1, 10):
            min_ovlp = i * 0.1

            # fisher test
            fisher = bed.fisher(bed_overlap_with, f=min_ovlp, genome=assembly)
            overlaps_count = fisher.table['in -a']['in -b']
            left_tail_fisher_pvalue = fisher.left_tail
            right_tail_fisher_pvalue = fisher.right_tail
            two_tail_fisher_pvalue = fisher.two_tail

            # jaccard index
            jaccard = bed.jaccard(bed_overlap_with, f=min_ovlp)
            jaccard_index = jaccard['jaccard']

            row_array = [bed_encyclopedia, biosample_name, bed_overlap_with_encyclopedia, a_size, b_size,
                         min_ovlp, overlaps_count, left_tail_fisher_pvalue, right_tail_fisher_pvalue,
                         two_tail_fisher_pvalue, jaccard_index]

            temp_df = pd.DataFrame([row_array], columns=columns)
            tests_df = tests_df.append(temp_df)

        tests_df.reset_index(inplace=True, drop=True)
        return tests_df
