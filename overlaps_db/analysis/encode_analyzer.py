from overlaps_db.analysis.overlap_analyzer import OverlapAnalyzer
from overlaps_db.store.hdf_store_manager import HdfStoreManager

import overlaps_db.utils.data_process_utils as utils
from overlaps_db.utils.timer import Timer


class EncodeAnalyzer(OverlapAnalyzer):
    def __init__(self, storage_path):
        super(EncodeAnalyzer, self).__init__(storage_path)

    def perform_overlap_analysis_with_fantom_permissive(self, assembly='hg19', method=None, overlap_intervals=10,
                                                        samples_num=20, biosample_type=None):
        storage_layer = HdfStoreManager(self.storage_path)

        fantom_file_name = utils.build_permissive_file_name()
        fantom_file_name_bed = utils.build_bed_file_name(fantom_file_name)
        fantom_bed = storage_layer.read_bed_file('fantom_staging.hdf', fantom_file_name_bed)
        fantom_bed_sorted = fantom_bed.sort()

        random_fantom_bed = self.build_random_from_bed(fantom_bed, assembly)
        shuffled_fantom_bed = self.build_shuffled_from_bed(fantom_bed, assembly)

        annotation_df = storage_layer.read_dataframe('encode_staging.hdf', 'encode_metadata')
        annotation_df.reset_index(inplace=True, drop=True)
        if assembly:
            annotation_df = annotation_df.query('assembly==@assembly')
        if method:
            annotation_df = annotation_df.query('method==@method')
        if biosample_type:
            annotation_df = annotation_df.query('biosample_type==@biosample_type')

        print(len(annotation_df), "experiments found")

        print("Computing fisher test, jaccard index for ENCODE and FANTOM overlaps,"
              " and z tests over random and shuffled null models...")
        tests_df = self.init_fisher_jaccard_z_tests_df()

        for index, row in annotation_df.iterrows():
            biosample_term_id = row['biosample_term_id']
            biosample_term_name = row['biosample_term_name']

            print("Processing", biosample_term_name)
            timer = Timer()
            timer.start()

            biosample_bed_name = utils.build_biosample_bed_file_name(biosample_term_id, assembly, method)
            biosample_bed = storage_layer.read_bed_file('encode_staging.hdf', biosample_bed_name)
            biosample_bed_sorted = biosample_bed.sort()

            tests_df = tests_df.append(
                self.compute_fisher_jaccard_z_tests(biosample_bed_sorted, fantom_bed_sorted, 'ENCODE', 'FANTOM',
                                                    biosample_term_name,
                                                    assembly, overlap_intervals, samples_num))

            tests_random_df = self.compute_fisher_jaccard_tests(biosample_bed_sorted, random_fantom_bed, 'ENCODE',
                                                                'RANDOM', biosample_term_name, assembly,
                                                                overlap_intervals)

            tests_shuffled_df = self.compute_fisher_jaccard_tests(biosample_bed_sorted, shuffled_fantom_bed, 'ENCODE',
                                                                  'SHUFFLED', biosample_term_name, assembly,
                                                                  overlap_intervals)

            tests_df = tests_df.append(tests_random_df)
            tests_df = tests_df.append(tests_shuffled_df)

            print(timer.elapsed())

        print("Storing in stats...")
        tests_df.reset_index(inplace=True, drop=True)
        storage_layer.store_dataframe(tests_df, 'stats.hdf', 'encode_fantom_tests')

        print("Completed")

    def perform_reldist_analsys_with_fantom_permissive(self, assembly='hg19', method=None, biosample_type=None):

        storage_layer = HdfStoreManager(self.storage_path)

        fantom_file_name = utils.build_permissive_file_name()
        fantom_file_name_bed = utils.build_bed_file_name(fantom_file_name)
        fantom_bed = storage_layer.read_bed_file('fantom_staging.hdf', fantom_file_name_bed)
        fantom_bed_sorted = fantom_bed.sort()

        random_fantom_bed = self.build_random_from_bed(fantom_bed, assembly)
        shuffled_fantom_bed = self.build_shuffled_from_bed(fantom_bed, assembly)

        annotation_df = storage_layer.read_dataframe('encode_staging.hdf', 'encode_metadata')
        annotation_df.reset_index(inplace=True, drop=True)
        if assembly:
            annotation_df = annotation_df.query('assembly==@assembly')
        if method:
            annotation_df = annotation_df.query('method==@method')
        if biosample_type:
            annotation_df = annotation_df.query('biosample_type==@biosample_type')

        print(len(annotation_df), "experiments found")

        print("Computing relative distance analysis for ENCODE and FANTOM overlaps")
        reldist_df = self.init_reldist_df()

        for index, row in annotation_df.iterrows():
            biosample_term_id = row['biosample_term_id']
            biosample_term_name = row['biosample_term_name']

            print("Processing", biosample_term_name)
            timer = Timer()
            timer.start()

            biosample_bed_name = utils.build_biosample_bed_file_name(biosample_term_id, assembly, method)
            biosample_bed = storage_layer.read_bed_file('encode_staging.hdf', biosample_bed_name)
            biosample_bed_sorted = biosample_bed.sort()

            reldist_df = reldist_df.append(
                self.compute_reldist(biosample_bed_sorted, fantom_bed_sorted, 'ENCODE', 'FANTOM',
                                     biosample_term_name)
            )

            reldist_df = reldist_df.append(
                self.compute_reldist(biosample_bed_sorted, random_fantom_bed, 'ENCODE', 'RANDOM',
                                     biosample_term_name)
            )

            reldist_df = reldist_df.append(
                self.compute_reldist(biosample_bed_sorted, shuffled_fantom_bed, 'ENCODE', 'SHUFFLED',
                                     biosample_term_name)
            )

            print(timer.elapsed())

        reldist_df.reset_index(inplace=True, drop=True)
        storage_layer.store_dataframe(reldist_df, 'stats.hdf', 'encode_fantom_reldist')
