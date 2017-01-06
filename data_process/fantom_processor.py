from data_process.processor import Processor
import pandas as pd
import os
from pybedtools import BedTool


class FantomBedProcessor(Processor):
    def __init__(self, download_path, staging_path):
        super(FantomBedProcessor, self).__init__(download_path, staging_path)

    def process_permissive(self, permissive_filename):
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

        output_filename = self.staging_path + "/permissive/"
        if not os.path.exists(output_filename):
            os.makedirs(output_filename)
        output_filename = output_filename + experiment + ".csv"

        bed_df.to_csv(output_filename, index=None, sep='\t')

        print("Completed")


# TODO test only
fantom_download_path = "/Users/manuel/development/thesis/download/FANTOM"
fantom_staging_path = "/Users/manuel/development/thesis/staging/FANTOM"
perm_filename = "permissive_enhancers.bed"

processor = FantomBedProcessor(fantom_download_path, fantom_staging_path)
processor.process_permissive(perm_filename)
