#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys
import pandas as pd
from overlaps_db.cli import cli
from pybedtools import BedTool

# TODO use this paths
# "/Users/manuel/development/thesis/download"
# "/Users/manuel/development/thesis/staging"
# "/Users/manuel/development/thesis/overlap"

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


class EncodeOverlapper:
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

    @staticmethod
    def build_filtered_file_name(assembly, method):
        encode_file_name = "filtered_"
        if assembly:
            encode_file_name += assembly
        if method:
            encode_file_name += method
        return encode_file_name

    def overlap_bed_files(self, bed_filepath, bed_overlap_with_filepath, min_overlap):
        bed = BedTool(bed_filepath)
        intesect_bed = BedTool(bed_overlap_with_filepath)

        bed_intersect_bed = bed.intersect(intesect_bed, loj=True, f=min_overlap)
        return bed_intersect_bed

    def overlap_filtered_with_fantom_permissive(self, assembly='hg19', method=None, min_overlap=0.1, export_temp=False):
        encode_file_name = self.build_filtered_file_name(assembly, method)

        encode_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".csv"
        # fantom_file_path = self.staging_path + "/FANTOM/permissive/PERMISSIVE.csv"
        encode_bed_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".bed"
        fantom_bed_filepath = self.staging_path + "/FANTOM/permissive/PERMISSIVE.bed"

        encode_bed = BedTool(encode_bed_file_path)
        fantom_bed = BedTool(fantom_bed_filepath)

        print("ENCODE bed file:", encode_bed_file_path)
        print("FANTOM bed file:", fantom_bed_filepath)
        print("Starting overlap...")

        encode_intersect_fantom = encode_bed.intersect(fantom_bed, loj=True, f=min_overlap)
        print(len(encode_intersect_fantom), "ENCODE.intersect(FANTOM) results")

        overlap_column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'FA_chrom', 'FA_start', 'FA_end',
                                'FA_name', 'FA_score', 'FA_strand']

        overlap_df = encode_intersect_fantom.to_dataframe(names=overlap_column_names)

        overlap_df['size'] = overlap_df.apply(lambda row: self.compute_size(row), axis=1)

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_FANTOM_overlapped_temp.csv"
            print("Exporting temporary overlapped file (no ENCODE details merge) to:", output_filename_temp)
            overlap_df.to_csv(output_filename_temp, index=None, sep='\t')

        print("Adding details from FANTOM...")

        overlap_df['FA_size'] = overlap_df.apply(lambda row: self.compute_size(row, 'FA'), axis=1)
        overlap_df['FA_method'] = overlap_df.apply(lambda row: self.set_col_value(row, 'FA_name', 'CAGE_TCs'), axis=1)
        overlap_df['FA_ovlp_len'] = overlap_df.apply(lambda row: self.compute_ovlp_len(row, 'FA'), axis=1)
        overlap_df['FA_ovlp_pct'] = overlap_df.apply(lambda row: self.compute_ovlp_pct(row, 'FA'), axis=1)
        overlap_df['FA_encyclopedia'] = overlap_df.apply(lambda row: self.set_col_value(row, 'FA_name', 'FANTOM'), axis=1)

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_FANTOM_overlapped_temp_detail.csv"
            print("Exporting temporary overlapped detailed file (no ENCODE details merge) to:", output_filename_temp)
            overlap_df.to_csv(output_filename_temp, index=None, sep='\t')

        print("ENCODE details file:", encode_file_path)
        encode_details_df = pd.read_csv(encode_file_path, sep='\t')

        print("Merging details from ENCODE...")
        overlap_full_detail_encode_df = overlap_df.merge(
            encode_details_df[['candidate_id', 'assembly', 'biosample_term_id', 'biosample_term_name', 'biosample_type',
                               'description', 'developmental_slims', 'encyclopedia',
                               'organ_slims', 'system_slims', 'method']],
            how='left', left_on='name', right_on='candidate_id'
        )

        print("Rearranging columns...")
        overlap_full_detail_encode_df = overlap_full_detail_encode_df[
            ['chrom', 'start', 'end', 'name', 'score', 'strand', 'size', 'method', 'description', 'assembly',
             'biosample_type',
             'biosample_term_id', 'biosample_term_name', 'developmental_slims', 'system_slims', 'organ_slims',
             'encyclopedia', 'FA_chrom', 'FA_start', 'FA_end',
             'FA_name', 'FA_score', 'FA_size', 'FA_method', 'FA_ovlp_len',
             'FA_ovlp_pct', 'FA_encyclopedia']]

        output_filename = self.overlap_path + "/" + encode_file_name + "_FANTOM_overlapped.csv"
        print("Exporting overlapped file to:", output_filename)
        overlap_full_detail_encode_df.to_csv(output_filename, index=None, sep='\t')

        print("Completed")

    def overlap_filtered_with_dbsuper(self, assembly='hg19', method=None, min_overlap=0.1, export_temp=False):
        encode_file_name = self.build_filtered_file_name(assembly, method)

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

        overlap_column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'SE_chrom', 'SE_start', 'SE_end',
                                'SE_name', 'SE_score']

        overlap_df = encode_intersect_dbsuper.to_dataframe(names=overlap_column_names)
        overlap_df['size'] = overlap_df.apply(lambda row: self.compute_size(row), axis=1)

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_dbSUPER_overlapped_temp.csv"
            print("Exporting temporary overlapped file (no ENCODE details merge) to:", output_filename_temp)
            overlap_df.to_csv(output_filename_temp, index=None, sep='\t')

        print("dbSUPER details file:", dbsuper_file_path)
        full_dbsuper_details_df = pd.read_csv(dbsuper_file_path, sep='\t', index_col='index')

        print("Merging details from dbSUPER...")
        overlap_full_detail_df = overlap_df.merge(
            full_dbsuper_details_df[['ID', 'Size', 'Associated Gene', 'Method', 'Rank', 'Cell/Tissue']],
            how='left', left_on='SE_name', right_on='ID')

        overlap_full_detail_df = overlap_full_detail_df.drop('ID', axis=1)
        overlap_full_detail_df = overlap_full_detail_df.drop('Rank', axis=1)

        overlap_full_detail_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'SE_chrom',
                                          'SE_start', 'SE_end', 'SE_name', 'SE_score', 'size', 'SE_size',
                                          'SE_associated_gene', 'SE_method', 'SE_biosample']

        overlap_full_detail_df['SE_ovlp_len'] = overlap_full_detail_df.apply(
            lambda row: self.compute_ovlp_len(row, 'SE'),
            axis=1)
        overlap_full_detail_df['SE_ovlp_pct'] = overlap_full_detail_df.apply(
            lambda row: self.compute_ovlp_pct(row, 'SE'),
            axis=1)

        overlap_full_detail_df['SE_associated_gene'].fillna('.', inplace=True)
        overlap_full_detail_df['SE_method'].fillna('.', inplace=True)
        overlap_full_detail_df['SE_biosample'].fillna('.', inplace=True)
        overlap_full_detail_df['SE_encyclopedia'] = 'dbSUPER'

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_dbSUPER_overlapped_temp_detail.csv"
            print("Exporting temporary overlapped detailed file (no ENCODE details merge) to:", output_filename_temp)
            overlap_full_detail_df.to_csv(output_filename_temp, index=None, sep='\t')

        print("ENCODE details file:", encode_file_path)
        encode_details_df = pd.read_csv(encode_file_path, sep='\t')

        print("Merging details from ENCODE...")
        overlap_full_detail_encode_df = overlap_full_detail_df.merge(
            encode_details_df[['candidate_id', 'assembly', 'biosample_term_id', 'biosample_term_name', 'biosample_type',
                               'description', 'developmental_slims', 'encyclopedia',
                               'organ_slims', 'system_slims', 'method']],
            how='left', left_on='name', right_on='candidate_id'
        )

        print("Rearranging columns...")
        overlap_full_detail_encode_df = overlap_full_detail_encode_df[
            ['chrom', 'start', 'end', 'name', 'score', 'strand', 'size', 'method', 'description', 'assembly',
             'biosample_type',
             'biosample_term_id', 'biosample_term_name', 'developmental_slims', 'system_slims', 'organ_slims',
             'encyclopedia', 'SE_chrom', 'SE_start', 'SE_end',
             'SE_name', 'SE_score', 'SE_size', 'SE_associated_gene', 'SE_method', 'SE_biosample', 'SE_ovlp_len',
             'SE_ovlp_pct', 'SE_encyclopedia']]

        output_filename = self.overlap_path + "/" + encode_file_name + "_dbSUPER_overlapped.csv"
        print("Exporting overlapped file to:", output_filename)
        overlap_full_detail_encode_df.to_csv(output_filename, index=None, sep='\t')

        print("Completed")
