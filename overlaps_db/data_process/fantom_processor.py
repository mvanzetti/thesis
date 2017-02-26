import os

from pybedtools import BedTool

from overlaps_db.data_process.processor import Processor
from overlaps_db.store.hdf_store_manager import HdfStoreManager


class FantomBedProcessor(Processor):
    def __init__(self, download_path, staging_path, storage_path):
        subfolder = "/FANTOM"
        super(FantomBedProcessor, self).__init__(download_path + subfolder, staging_path + subfolder, storage_path)

    def process_permissive(self):
        storage_layer = HdfStoreManager(self.storage_path)

        permissive_filename = "permissive_enhancers.bed"
        bed_filepath = self.download_path + "/" + permissive_filename
        bed_file = BedTool(bed_filepath)
        bed_df = bed_file.to_dataframe()

        print("Found", len(bed_df), "permissive enhancers")

        print("Processing and adding info...")
        bed_df['seq_id'] = range(bed_df.index.size)
        bed_df['seq_id'] = bed_df['seq_id'].astype(str)

        assembly = 'hg19'
        encyclopedia = 'FANTOM'
        version = 5
        description = 'Permissive enhancers from http://enhancer.binf.ku.dk/presets/'
        organism = '/organisms/human'
        experiment = 'PERMISSIVE'
        method = 'CAGE_TCs'

        print("Assigning the following defaults")
        print("     assembly:", assembly)
        print("     encyclopedia:", encyclopedia)
        print("     version:", version)
        print("     description:", description)
        print("     organism:", organism)
        print("     experiment:", experiment)
        print("     method:", method)

        bed_df['candidate_id'] = encyclopedia + '.' + str(version) + '.' + experiment + '.' + bed_df['seq_id']

        bed_df['assembly'] = assembly
        bed_df['biosample_term_id'] = ''
        bed_df['biosample_term_name'] = ''
        bed_df['biosample_type'] = ''
        bed_df['description'] = description
        bed_df['developmental_slims'] = ''
        bed_df['encyclopedia'] = encyclopedia
        bed_df['encyclopedia_version'] = version
        bed_df['month_released'] = ''
        bed_df['organ_slims'] = ''
        bed_df['organism'] = organism
        bed_df['system_slims'] = ''
        bed_df['method'] = method

        bed_df = bed_df.drop('seq_id', axis=1)
        bed_df = bed_df.drop('itemRgb', axis=1)

        # output_filename_path = self.staging_path + "/permissive/"
        # if not os.path.exists(output_filename_path):
        #     os.makedirs(output_filename_path)
        # output_filename = output_filename_path + experiment + ".csv"
        # output_bed_filename = output_filename_path + experiment + ".bed"

        output_filename = 'permissive'
        print("Exporting permissive file to staging:", output_filename)
        storage_layer.store_dataframe_fixed(bed_df, 'fantom_staging.hdf', output_filename)

        # print("Exporting permissive file to:", output_filename)
        # bed_df.to_csv(output_filename, index=None, sep='\t')

        output_bed_filename = output_filename + '_bed'
        print("Exporting permissive file in bed format (substituting names with candidate_ids) to staging:",
              output_bed_filename)
        bed_df['name'] = bed_df['candidate_id']
        bed_df = bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]

        storage_layer.store_dataframe_fixed(bed_df, 'fantom_staging.hdf', output_bed_filename)
        # bed_df.to_csv(output_bed_filename, index=None, sep='\t', header=None)

        print("Completed")
