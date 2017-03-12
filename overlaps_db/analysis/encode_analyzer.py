from overlaps_db.analysis.overlap_analyzer import OverlapAnalyzer
from overlaps_db.store.hdf_store_manager import HdfStoreManager

import overlaps_db.utils.data_process_utils as utils
from overlaps_db.utils.timer import Timer


class EncodeAnalyzer(OverlapAnalyzer):
    def __init__(self, storage_path):
        super(EncodeAnalyzer, self).__init__(storage_path)

    def build_table_name(self, name, assembly, method, biosample_type, sample_param):
        table_name = name

        if assembly:
            table_name += "_" + assembly
        if method:
            table_name += "_" + method
        if biosample_type:
            table_name += "_" + biosample_type.replace(" ", "_")
        if sample_param:
            table_name += "_" + sample_param

        return table_name

    def read_annotation_df(self, assembly, method, biosample_type):
        storage_layer = HdfStoreManager(self.storage_path)
        annotation_df = storage_layer.read_dataframe('encode_staging.hdf', 'encode_metadata')
        annotation_df.drop_duplicates(inplace=True)
        annotation_df.reset_index(inplace=True, drop=True)

        if assembly:
            annotation_df = annotation_df.query('assembly==@assembly')
        if method:
            annotation_df = annotation_df.query('method==@method')
        if biosample_type:
            annotation_df = annotation_df.query('biosample_type==@biosample_type')

        return annotation_df

    def read_biosample_bed(self, biosample_bed_name):
        storage_layer = HdfStoreManager(self.storage_path)
        return storage_layer.read_bed_file('encode_staging.hdf', biosample_bed_name)

    def perform_overlap_analysis_with(self, sample_bed_with, sample_bed_with_name, assembly, method,
                                      overlap_intervals,
                                      samples_num, biosample_type, avoid_z_test=False):

        sample_bed_sorted = sample_bed_with.sort()
        random_bed = self.build_random_from_bed(sample_bed_with, assembly)
        shuffled_bed = self.build_shuffled_from_bed(sample_bed_with, assembly)

        table_name = self.build_table_name('overlap', assembly, method, biosample_type, sample_bed_with_name)
        annotation_df = self.read_annotation_df(assembly, method, biosample_type)
        print(len(annotation_df), "experiments found")

        print("Computing fisher test, jaccard index for ENCODE and", sample_bed_with_name,
              "overlaps and z tests over random and shuffled null models...")
        tests_df = self.init_fisher_jaccard_z_tests_df()

        for index, row in annotation_df.iterrows():
            bio_type = row['biosample_type']
            biosample_term_id = row['biosample_term_id']
            biosample_term_name = row['biosample_term_name']

            print("Processing", biosample_term_name)
            timer = Timer()
            timer.start()

            biosample_bed_name = utils.build_biosample_bed_file_name(biosample_term_id, assembly, method)
            biosample_bed = self.read_biosample_bed(biosample_bed_name)
            biosample_bed_sorted = biosample_bed.sort()

            tests_df = tests_df.append(
                self.compute_fisher_jaccard_z_tests(biosample_bed_sorted, sample_bed_sorted, 'ENCODE',
                                                    sample_bed_with_name,
                                                    bio_type, biosample_term_name,
                                                    assembly, overlap_intervals, samples_num, avoid_z_test))

            tests_random_df = self.compute_fisher_jaccard_tests(biosample_bed_sorted, random_bed, 'ENCODE',
                                                                'RANDOM', bio_type, biosample_term_name, assembly,
                                                                overlap_intervals)

            tests_shuffled_df = self.compute_fisher_jaccard_tests(biosample_bed_sorted, shuffled_bed, 'ENCODE',
                                                                  'SHUFFLED', bio_type, biosample_term_name, assembly,
                                                                  overlap_intervals)

            tests_df = tests_df.append(tests_random_df)
            tests_df = tests_df.append(tests_shuffled_df)

            print(timer.elapsed())

        return tests_df, table_name

    def perform_reldist_analsys_with(self, sample_bed_with, sample_bed_with_name, assembly, method, biosample_type):

        sample_bed_sorted = sample_bed_with.sort()

        random_bed = self.build_random_from_bed(sample_bed_with, assembly)
        shuffled_bed = self.build_shuffled_from_bed(sample_bed_with, assembly)

        annotation_df = self.read_annotation_df(assembly, method, biosample_type)
        table_name = self.build_table_name('reldist', assembly, method, biosample_type, sample_bed_with_name)

        print(len(annotation_df), "experiments found")

        print("Computing relative distance analysis for ENCODE and", sample_bed_with_name, "overlaps")
        reldist_df = self.init_reldist_df()

        for index, row in annotation_df.iterrows():
            bio_type = row['biosample_type']
            biosample_term_id = row['biosample_term_id']
            biosample_term_name = row['biosample_term_name']

            print("Processing", biosample_term_name)
            timer = Timer()
            timer.start()

            biosample_bed_name = utils.build_biosample_bed_file_name(biosample_term_id, assembly, method)
            biosample_bed = self.read_biosample_bed(biosample_bed_name)
            biosample_bed_sorted = biosample_bed.sort()

            try:
                reldist_df = reldist_df.append(
                    self.compute_reldist(biosample_bed_sorted, sample_bed_sorted, 'ENCODE', sample_bed_with_name,
                                         bio_type, biosample_term_name)
                )

                reldist_df = reldist_df.append(
                    self.compute_reldist(biosample_bed_sorted, random_bed, 'ENCODE', 'RANDOM', bio_type,
                                         biosample_term_name)
                )

                reldist_df = reldist_df.append(
                    self.compute_reldist(biosample_bed_sorted, shuffled_bed, 'ENCODE', 'SHUFFLED', bio_type,
                                         biosample_term_name)
                )
            except:
                print("Error")

            print(timer.elapsed())

        return reldist_df, table_name

    # def perform_overlap_analysis_with_repeatmasker(self, assembly='hg19', method=None, overlap_intervals=10,
    #                                                samples_num=20, biosample_type=None, repeat_class=None):
    #
    #     storage_layer = HdfStoreManager(self.storage_path)
    #
    #     repeat_class_file_name = utils.build_repeatmasker_file_name(repeat_class)
    #     repeat_class_file_name_bed = utils.build_bed_file_name(repeat_class_file_name)
    #     repeat_bed = storage_layer.read_bed_file('repeatmasker_staging.hdf', repeat_class_file_name_bed)
    #
    #     tests_df, table_name = self.perform_overlap_analysis_with(repeat_bed, repeat_class_file_name, assembly, method,
    #                                                               overlap_intervals, samples_num, biosample_type,
    #                                                               avoid_z_test=True)
    #     print("Storing in stats...")
    #     tests_df.reset_index(inplace=True, drop=True)
    #
    #     storage_layer.store_dataframe(tests_df, 'encode_repeatmasker_stats.hdf', table_name,
    #                                   queryable_cols=['biosample_name'])
    #
    #     print("Completed")

    def perform_overlap_analysis_with_fantom_permissive(self, assembly='hg19', method=None, overlap_intervals=10,
                                                        samples_num=20, biosample_type=None):
        storage_layer = HdfStoreManager(self.storage_path)

        fantom_file_name = utils.build_permissive_file_name()
        fantom_file_name_bed = utils.build_bed_file_name(fantom_file_name)
        fantom_bed = storage_layer.read_bed_file('fantom_staging.hdf', fantom_file_name_bed)

        tests_df, table_name = self.perform_overlap_analysis_with(fantom_bed, fantom_file_name, assembly, method,
                                                                  overlap_intervals, samples_num, biosample_type)
        print("Storing in stats...")
        tests_df.reset_index(inplace=True, drop=True)

        storage_layer.store_dataframe(tests_df, 'encode_fantom_stats.hdf', table_name,
                                      queryable_cols=['biosample_name'])

        print("Completed")

    # def perform_reldist_analsys_with_repeatmasker(self, assembly='hg19', method=None, biosample_type=None,
    #                                               repeat_class=None):
    #
    #     storage_layer = HdfStoreManager(self.storage_path)
    #
    #     repeat_class_file_name = utils.build_repeatmasker_file_name(repeat_class)
    #     repeat_class_file_name_bed = utils.build_bed_file_name(repeat_class_file_name)
    #     repeat_bed = storage_layer.read_bed_file('repeatmasker_staging.hdf', repeat_class_file_name_bed)
    #
    #     reldist_df, table_name = self.perform_reldist_analsys_with(repeat_bed, repeat_class_file_name, assembly, method,
    #                                                                biosample_type)
    #
    #     print("Storing in stats...")
    #     reldist_df.reset_index(inplace=True, drop=True)
    #
    #     storage_layer.store_dataframe(reldist_df, 'encode_repeatmasker_stats.hdf', table_name,
    #                                   queryable_cols=['biosample_name'])
    #
    #     print("Completed")

    def perform_reldist_analsys_with_fantom_permissive(self, assembly='hg19', method=None, biosample_type=None):

        storage_layer = HdfStoreManager(self.storage_path)

        fantom_file_name = utils.build_permissive_file_name()
        fantom_file_name_bed = utils.build_bed_file_name(fantom_file_name)
        fantom_bed = storage_layer.read_bed_file('fantom_staging.hdf', fantom_file_name_bed)

        reldist_df, table_name = self.perform_reldist_analsys_with(fantom_bed, fantom_file_name, assembly, method,
                                                                   biosample_type)

        print("Storing in stats...")
        reldist_df.reset_index(inplace=True, drop=True)
        storage_layer.store_dataframe(reldist_df, 'encode_fantom_stats.hdf', table_name,
                                      queryable_cols=['biosample_name'])
        print("Completed")
