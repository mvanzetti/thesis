import os
import sys

import pandas as pd
from pybedtools import BedTool

from overlaps_db.data_process.processor import Processor
from overlaps_db.store.hdf_store_manager import HdfStoreManager

import utils.data_process_utils as utils


class EncodeBedProcessor(Processor):
    def __init__(self, download_path, staging_path, storage_path):
        subfolder = "/ENCODE"
        super(EncodeBedProcessor, self).__init__(download_path + subfolder, staging_path + subfolder, storage_path)

    @staticmethod
    def label_method(row):
        description = row['description']
        if 'DNase-only' in description:
            return 'DNase'
        if 'H3K27ac-only' in description:
            return 'H3K27ac'
        if 'DNase and H3K27ac' in description:
            return 'DNase_H3K27ac'
        else:
            return 'Unknown'

    def prepare(self):
        storage_layer = HdfStoreManager(self.storage_path)

        print("Importing main annotation file...")
        # file_list_name = self.download_path + "/" + annotation_filename
        # file_list_output_name = self.staging_path + "/" + prepared_filename
        # annotation_df = pd.read_csv(file_list_name, sep='\t', index_col='index')
        annotation_df = storage_layer.read_dataframe('downloads.hdf', 'encode_metadata')

        imported_df = annotation_df.query('imported == True')
        count_downloaded = len(imported_df)
        count_total = len(annotation_df)
        print("Focusing on the", count_downloaded, "imported entries out of", count_total)

        print("Building .bed.gz files references...")
        output_df = imported_df[[
            'accession',
            'assembly',
            'biosample_term_id',
            'biosample_term_name',
            'biosample_type',
            'description',
            'developmental_slims',
            'encyclopedia_version',
            'month_released',
            'organ_slims',
            'organism',
            'system_slims'
        ]]

        output_df.insert(output_df.columns.size, 'bed_filename', "")
        output_df.insert(output_df.columns.size, 'bed_filepath', "")
        output_df.insert(output_df.columns.size, 'merged', False)
        output_df.insert(output_df.columns.size, 'method', output_df.apply(lambda r: self.label_method(r), axis=1))

        # output_df['bed_filename'] = ""
        # output_df['bed_filepath'] = ""
        # output_df['merged'] = False
        # output_df['method'] = output_df.apply(lambda r: self.label_method(r), axis=1)

        # build a df with info and path to bed.gz file
        for index, row in output_df.iterrows():
            accession = output_df.get_value(index, 'accession')
            folder_path = self.download_path + "/" + accession + "/files"
            for dir in os.listdir(folder_path):
                current_dir = folder_path + "/" + dir

                if not os.path.exists(current_dir + "/"):
                    continue

                files = os.listdir(current_dir)
                for file in files:
                    if file.endswith(".bed.gz"):
                        output_df.set_value(index, 'bed_filename', file)
                        file_path = current_dir + "/" + file
                        output_df.set_value(index, 'bed_filepath', file_path)

        print("Exporting metadata and .bed.gz file references to staging")
        utils.stringify_columns(output_df, ['developmental_slims', 'organ_slims', 'system_slims'])
        storage_layer.store_dataframe_fixed(output_df, 'encode_staging.hdf', 'encode_metadata')
        # output_df.to_csv(file_list_output_name, sep='\t')
        print("Completed")

    def process_bed_file(self, index, annotation_df):
        bed_filepath = annotation_df.get_value(index, 'bed_filepath')
        bed_file = BedTool(bed_filepath)
        bed_df = bed_file.to_dataframe()

        # assembly_val = annotation_df.get_value(index, 'assembly')
        # assembly_list = ast.literal_eval(assembly_val)
        # assembly = assembly_list[0] if len(assembly_list) > 0 else ''

        bed_df['seq_id'] = range(bed_df.index.size)
        bed_df['seq_id'] = bed_df['seq_id'].astype(str)

        encyclopedia = 'ENCODE'
        version = annotation_df.get_value(index, 'encyclopedia_version')
        experiment = annotation_df.get_value(index, 'bed_filename').split(sep='.')[0]

        bed_df['candidate_id'] = encyclopedia + '.' + str(version) + '.' + experiment + '.' + bed_df['seq_id']

        bed_df['assembly'] = annotation_df.get_value(index, 'assembly')
        bed_df['biosample_term_id'] = annotation_df.get_value(index, 'biosample_term_id')
        bed_df['biosample_term_name'] = annotation_df.get_value(index, 'biosample_term_name')
        bed_df['biosample_type'] = annotation_df.get_value(index, 'biosample_type')
        bed_df['description'] = annotation_df.get_value(index, 'description')
        bed_df['developmental_slims'] = annotation_df.get_value(index, 'developmental_slims')
        bed_df['encyclopedia'] = encyclopedia
        bed_df['encyclopedia_version'] = version
        bed_df['month_released'] = annotation_df.get_value(index, 'month_released')
        bed_df['organ_slims'] = annotation_df.get_value(index, 'organ_slims')
        bed_df['organism'] = annotation_df.get_value(index, 'organism')
        bed_df['system_slims'] = annotation_df.get_value(index, 'system_slims')
        bed_df['method'] = annotation_df.get_value(index, 'method')

        bed_df = bed_df.drop('seq_id', axis=1)
        bed_df = bed_df.drop('thickStart', axis=1)
        bed_df = bed_df.drop('thickEnd', axis=1)
        bed_df = bed_df.drop('itemRgb', axis=1)

        # output_filename = self.staging_path + "/merged/" + experiment + ".csv"

        return bed_df

    def filter(self, assembly=None, method=None, force=False):
        storage_layer = HdfStoreManager(self.storage_path)
        # prepared_filename = "enhancer-like-bed_refs.csv"
        # file_list_name = self.staging_path + "/" + prepared_filename
        # annotation_df = pd.read_csv(file_list_name, sep='\t', index_col='index')
        annotation_df = storage_layer.read_dataframe('encode_staging.hdf', 'encode_metadata')
        annotation_df.reset_index(inplace=True, drop=True)

        filter_file_name = utils.build_filtered_file_name(assembly, method)

        if assembly:
            annotation_df = annotation_df.query('assembly==@assembly')
        if method:
            annotation_df = annotation_df.query('method==@method')

        print(len(annotation_df), "experiments found")
        print("Building filtered file in staging...")

        full_bed_df = pd.DataFrame()

        for index, row in annotation_df.iterrows():
            try:
                if annotation_df.get_value(index, 'merged') and not force:
                    continue

                if pd.isnull(annotation_df.get_value(index, 'bed_filepath')):
                    print("Accession", annotation_df.get_value(index, "accession"), "at index", index,
                          "has no bed file: skip")
                    continue

                bed_df = self.process_bed_file(index, annotation_df)
                biosample_bed_filename = utils.build_bed_file_name(
                    utils.build_biosample_file_name(row['biosample_term_id'], assembly, method))

                print("Exporting single biosample filtered file in bed format"
                      " (substituting names with candidate_ids) to staging:", biosample_bed_filename)

                biosample_bed_df = bed_df.copy()
                biosample_bed_df['name'] = biosample_bed_df['candidate_id']
                biosample_bed_df = biosample_bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
                storage_layer.store_dataframe_fixed(biosample_bed_df, 'encode_staging.hdf', biosample_bed_filename)

                full_bed_df = full_bed_df.append(bed_df)

            except ValueError as e:
                print("ValueError at index", index, ":{0}".format(e))
                print("Jumping to next row...")

            except:
                print("Unexpected error:", sys.exc_info()[0])
                print("Updating info file before to raise exception...")
                storage_layer.store_dataframe_fixed(annotation_df, 'encode_staging.hdf', 'encode_metadata')
                # annotation_df.to_csv(file_list_name, sep='\t')
                raise

        # output_filename = self.staging_path + "/filtered/" + filter_file_name + ".csv"
        # output_bed_filename = self.staging_path + "/filtered/" + filter_file_name + ".bed"

        print("Exporting filtered file to staging:", filter_file_name,
              "- queryable columns are [biosample_term_name, biosample_type]")
        storage_layer.store_dataframe(full_bed_df, 'encode_staging.hdf', filter_file_name,
                                      ['biosample_term_name', 'biosample_type'])
        # full_bed_df.to_csv(output_filename, index=None, sep='\t')

        output_bed_filename = filter_file_name + '_bed'
        print("Exporting filtered file in bed format (substituting names with candidate_ids) to staging:",
              output_bed_filename)
        full_bed_df['name'] = full_bed_df['candidate_id']
        full_bed_df_bed = full_bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]

        storage_layer.store_dataframe_fixed(full_bed_df_bed, 'encode_staging.hdf', output_bed_filename)

        # full_bed_df_bed.to_csv(output_bed_filename, index=None, sep='\t', header=None)

        print("Completed")
