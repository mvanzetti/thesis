import os
import sys

import pandas as pd
from pybedtools import BedTool

from overlaps_db.data_process.processor import Processor


class EpigenomicsRoadmapProcessor(Processor):
    def __init__(self, download_path, staging_path):
        subfolder = "/EpigenomicsRoadmap"
        super(EpigenomicsRoadmapProcessor, self).__init__(download_path + subfolder, staging_path + subfolder)

    @staticmethod
    def make_directory(directory):
        if not os.path.exists(directory):
            print("Directory", directory, "does not exist, creating...")
            os.makedirs(directory)

    def process_bed_file(self, index, metadata_df):
        eid = metadata_df.get_value(index, 'EID')
        bed_filepath = self.download_path + "/regions_enh_" + eid + ".bed.gz"

        if not os.path.exists(bed_filepath):
            print("Skipping: no bed file for", eid)
            return None

        bed_file = BedTool(bed_filepath)
        bed_df = bed_file.to_dataframe()

        bed_df['seq_id'] = range(bed_df.index.size)
        bed_df['seq_id'] = bed_df['seq_id'].astype(str)

        encyclopedia = 'ROADMAP'
        bed_df['candidate_id'] = encyclopedia + '.' + eid + '.' + bed_df['seq_id']

        bed_df = bed_df.drop('seq_id', axis=1)

        return bed_df

    def process(self):
        print("Copying metadata to staging...")
        read_metadata_filename = self.download_path + "/" + "roadmap_metadata.csv"
        metadata_df = pd.read_csv(read_metadata_filename, sep='\t')

        self.make_directory(self.staging_path)
        metadata_filename = self.staging_path + "/" + "roadmap_metadata.csv"
        metadata_df.to_csv(metadata_filename, sep='\t', index=None)

        processed_filename = "processed"

        print(len(metadata_df), "epigenomes found")
        print("Building unified bed file in staging...")

        full_bed_df = pd.DataFrame()

        for index, row in metadata_df.iterrows():
            bed_df = self.process_bed_file(index, metadata_df)
            if bed_df is None:
                continue
            full_bed_df = full_bed_df.append(bed_df)

        output_bed_filename = self.staging_path + "/" + processed_filename + ".bed"

        print("Exporting processed file in bed format (substituting names with candidate_ids) to:", output_bed_filename)
        full_bed_df['name'] = full_bed_df['candidate_id']
        full_bed_df_bed = full_bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
        full_bed_df_bed.to_csv(output_bed_filename, index=None, sep='\t', header=None)

        print("Completed")
