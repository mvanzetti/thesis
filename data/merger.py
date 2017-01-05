import sys
import pandas as pd
import os
import ast
from pybedtools import BedTool


class EncodeBedMerger:
    def __init__(self, download_path, staging_path):
        self.download_path = download_path
        self.staging_path = staging_path

    def prepare(self, annotation_filename, prepared_filename):
        print("Importing main annotation file...")
        file_list_name = self.download_path + "/" + annotation_filename
        file_list_output_name = self.staging_path + "/" + prepared_filename
        annotation_df = pd.read_csv(file_list_name, sep='\t', index_col='index')

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

        output_df['bed_filename'] = ""
        output_df['bed_filepath'] = ""
        output_df['merged'] = False

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

                        # bed_file = BedTool(file_path)
                        # bed_file

        print("Exporting .bed.gz file references...")
        output_df.to_csv(file_list_output_name, sep='\t')
        print("Completed")

    def merge(self, prepared_filename, force):
        file_list_name = self.staging_path + "/" + prepared_filename
        annotation_df = pd.read_csv(file_list_name, sep='\t', index_col='index')

        print("Merging annotations and info files in staging...")

        for index, row in annotation_df.iterrows():
            try:
                if annotation_df.get_value(index, 'merged') and not force:
                    continue

                if pd.isnull(annotation_df.get_value(index, 'bed_filepath')):
                    print("Accession", annotation_df.get_value(index, "accession"), "at index", index,
                          "has no bed file: skip")
                    continue

                bed_filepath = annotation_df.get_value(index, 'bed_filepath')
                bed_file = BedTool(bed_filepath)
                bed_df = bed_file.to_dataframe()

                assembly_val = annotation_df.get_value(index, 'assembly')
                assembly_list = ast.literal_eval(assembly_val)
                assembly = assembly_list[0] if len(assembly_list) > 0 else ''

                bed_df['seq_id'] = range(bed_df.index.size)
                bed_df['seq_id'] = bed_df['seq_id'].astype(str)

                encyclopedia = 'ENCODE'
                version = annotation_df.get_value(index, 'encyclopedia_version')
                experiment = annotation_df.get_value(index, 'bed_filename').split(sep='.')[0]

                bed_df['candidate_id'] = encyclopedia + '.' + str(version) + '.' + experiment + '.' + bed_df['seq_id']

                bed_df['assembly'] = assembly
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

                bed_df = bed_df.drop('seq_id', axis=1)
                bed_df = bed_df.drop('thickStart', axis=1)
                bed_df = bed_df.drop('thickEnd', axis=1)
                bed_df = bed_df.drop('itemRgb', axis=1)

                output_filename = self.staging_path + "/merged/" + experiment + ".csv"

                bed_df.to_csv(output_filename, index=None, sep='\t')
                annotation_df.set_value(index, 'merged', True)

            except ValueError as e:
                print("ValueError at index", index, ":{0}".format(e))
                print("Jumping to next row...")

            except:
                print("Unexpected error:", sys.exc_info()[0])
                print("Updating info file before to raise exception...")
                annotation_df.to_csv(file_list_name, sep='\t')
                raise

        print("Updating info file...")
        annotation_df.to_csv(file_list_name, sep='\t')
        print("Completed")

# # TODO test only
encode_path = "/Users/manuel/development/thesis/download/ENCODE"
encode_staging_path = "/Users/manuel/development/thesis/staging/ENCODE"
input_csv = "enhancer-like-annotations.csv"
output_csv = "enhancer-like-bed_refs.csv"

merger = EncodeBedMerger(encode_path, encode_staging_path)
# merger.prepare(input_csv, output_csv)
merger.merge(output_csv, force=False)
