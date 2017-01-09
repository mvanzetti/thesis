import pandas as pd

from pybedtools import BedTool


# TODO use this paths
# "/Users/manuel/development/thesis/download"
# "/Users/manuel/development/thesis/staging"
# "/Users/manuel/development/thesis/overlap"

class EncodeOverlapper:
    def __init__(self, download_path, staging_path, overlap_path):
        self.download_path = download_path
        self.staging_path = staging_path
        self.overlap_path = overlap_path

    def compute_ovlp_len(self, row):
        if row['SE_name'] == '.':
            return 0
        value = min(row['end'], row['SE_end']) - max(row['start'], row['SE_start'])
        return value

    def compute_ovlp_pct(self, row):
        size = abs(row['end'] - row['start'])
        value = row['SE_ovlp_len'] / size * 100.0
        return value

    def overlap_filtered_with_dbsuper(self, assembly='hg19', method=None, min_overlap=0.1):
        encode_file_name = "filtered_"
        if assembly:
            encode_file_name += assembly
        if method:
            encode_file_name += method

        encode_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".csv"
        dbsuper_file_path = self.download_path + "/dbSUPER/super-enhancers-annotations.csv"
        encode_bed_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".bed"
        dbsuper_bed_filepath = self.download_path + "/dbSUPER/" + "all_" + assembly + "_bed.bed"

        encode_bed = BedTool(encode_bed_file_path)
        dbsuper_bed = BedTool(dbsuper_bed_filepath)

        print("ENCODE bed file:", encode_bed_file_path)
        print("dbSUPER bed file:", dbsuper_bed_filepath)
        print("Starting overlap...")

        # full_dbSUPER_intersect_ENCODE = dbsuper_bed.intersect(encode_bed, wa=True, f=0.1)
        encode_intersect_dbsuper = encode_bed.intersect(dbsuper_bed, loj=True, f=min_overlap)

        # print(len(full_dbSUPER_intersect_ENCODE), "dbSUPER.intersect(ENCODE) results")
        print(len(encode_intersect_dbsuper), "ENCODE.intersect(dbSUPER) results")

        overlap_column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                                'SE_chrom', \
                                'SE_start', 'SE_end', 'SE_name', 'SE_score']
        overlap_df = encode_intersect_dbsuper.to_dataframe(names=overlap_column_names)

        print("dbSUPER details file:", dbsuper_file_path)
        full_dbsuper_details_df = pd.read_csv(dbsuper_file_path, sep='\t', index_col='index')

        print("Merging details from dbSUPER...")
        overlap_full_detail_df = overlap_df.merge(
            full_dbsuper_details_df[['ID', 'Size', 'Associated Gene', 'Method', 'Rank', 'Cell/Tissue']],
            how='left', left_on='SE_name', right_on='ID')

        overlap_full_detail_df = overlap_full_detail_df.drop('thickStart', axis=1)
        overlap_full_detail_df = overlap_full_detail_df.drop('thickEnd', axis=1)
        overlap_full_detail_df = overlap_full_detail_df.drop('itemRgb', axis=1)
        overlap_full_detail_df = overlap_full_detail_df.drop('ID', axis=1)
        overlap_full_detail_df = overlap_full_detail_df.drop('Rank', axis=1)

        overlap_full_detail_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'SE_chrom',
                                          'SE_start', 'SE_end', 'SE_name', 'SE_score', 'SE_size',
                                          'SE_associated_gene', 'SE_method', 'SE_biosample']

        overlap_full_detail_df['SE_ovlp_len'] = overlap_full_detail_df.apply(lambda row: self.compute_ovlp_len(row),
                                                                             axis=1)
        overlap_full_detail_df['SE_ovlp_pct'] = overlap_full_detail_df.apply(lambda row: self.compute_ovlp_pct(row),
                                                                             axis=1)

        overlap_full_detail_df['SE_associated_gene'].fillna('.', inplace=True)
        overlap_full_detail_df['SE_method'].fillna('.', inplace=True)
        overlap_full_detail_df['SE_biosample'].fillna('.', inplace=True)

        print("ENCODE details file:", encode_file_path)
        encode_details_df = pd.read_csv(encode_file_path, sep='\t')

        print("Merging details from ENCODE...")
        overlap_full_detail_encode_df = overlap_full_detail_df.merge(
            encode_details_df[['candidate_id', 'assembly', 'biosample_term_id', 'biosample_term_name', 'biosample_type',
                               'description', 'developmental_slims', 'encyclopedia', 'encyclopedia_version',
                               'organ_slims', 'system_slims', 'method']],
            how='left', left_on='name', right_on='candidate_id'
        )

        output_filename = self.overlap_path + "/" + encode_file_name + "_dbSUPER_overlapped.csv"
        print("Exporting overlapped file to:", output_filename)
        overlap_full_detail_encode_df.to_csv(output_filename, index=None, sep='\t')


# TODO test only

d = "/Users/manuel/development/thesis/download"
s = "/Users/manuel/development/thesis/staging"
o = "/Users/manuel/development/thesis/overlap"
overlapper = EncodeOverlapper(d, s, o)
overlapper.overlap_filtered_with_dbsuper(assembly='hg19', method='DNase_H3K27ac', min_overlap=0.1)
