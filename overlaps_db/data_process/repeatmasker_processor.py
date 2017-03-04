import os

from pybedtools import BedTool
import pandas as pd
from overlaps_db.data_process.processor import Processor
from overlaps_db.store.hdf_store_manager import HdfStoreManager
import overlaps_db.utils.data_process_utils as utils


class RepeatMaskerProcessor(Processor):
    def __init__(self, download_path, staging_path, storage_path):
        subfolder = "/RepeatMasker"
        super(RepeatMaskerProcessor, self).__init__(download_path + subfolder, staging_path + subfolder, storage_path)

    def process(self, assembly="hg19"):
        storage_layer = HdfStoreManager(self.storage_path)

        filename = assembly + ".fa.out"

        repeat_df = pd.read_csv(self.download_path + "/" + filename, sep=r"\s+", skiprows=2, header=None,
                                index_col=False)

        repeat_df.columns = ['SW_score', 'perc_div', 'perc_del', 'perc_ins', 'query_sequence',
                             'pos_in_query_begin', 'pos_in_query_end', 'pos_in_query_left', 'match_with_seq',
                             'matching_repeat', 'repeat_class_family',
                             'pos_in_repeat_begin', 'pos_in_repeat_end', 'pos_in_repeat_left', 'ID']

        nan_ID_indexes = repeat_df.isnull().query("ID == True").index.values

        if nan_ID_indexes.size > 0:
            print("Found", nan_ID_indexes.size, "nan IDs: filling...")
            for idx in nan_ID_indexes:
                new_id = max(repeat_df.ID) + 1
                repeat_df.set_value(idx, 'ID', new_id)
                print("Assigned ID", new_id, "to row at index", idx)

        repeat_df['ID'] = repeat_df.ID.astype(int).astype(str)

        print("Removing haplotype chromosomes, unplaced contigs and unlocalized contigs...")

        chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chrX',
                      'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                      'chr15', 'chr16', 'chr17', 'chr18', 'chr20', 'chrY', 'chr19',
                      'chr22', 'chr21']

        repeat_df = repeat_df.query("query_sequence in @chrom_list")

        repeat_classes = repeat_df.repeat_class_family.unique()

        print("Found", repeat_classes.size, "repeat classes")

        for class_name in repeat_classes:
            print("Processing", class_name, "...")
            class_df = repeat_df.query("repeat_class_family == @class_name")
            class_df.reset_index(inplace=True, drop=True)

            bed_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
            bed_df.chrom = class_df.query_sequence
            bed_df.start = class_df.pos_in_query_begin
            bed_df.end = class_df.pos_in_query_end
            bed_df.name = "RepeatMasker" + "." + pd.Series(
                class_df.index.values.astype(str)) + "." + class_df.matching_repeat + "." + class_df.ID
            bed_df.score = class_df.SW_score
            bed_df.strand = "."

            table_name = class_name.replace('/', '_').replace('?', '_qm').replace('-', '_')
            storage_layer.store_dataframe_fixed(class_df, 'repeatmasker_staging.hdf', table_name)

            bed = BedTool.from_dataframe(bed_df)
            bed_name = utils.build_bed_file_name(table_name)

            storage_layer.store_bed_file(bed, 'repeatmasker_staging.hdf', bed_name)

        print("Completed")
