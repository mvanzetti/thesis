#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys
import pandas as pd
from overlaps_db.cli import cli
from pybedtools import BedTool

from overlaps_db.data_process.overlapper import Overlapper
from overlaps_db.store.hdf_store_manager import HdfStoreManager
import overlaps_db.utils.data_process_utils as utils

# TODO use this paths
# "/Users/manuel/development/thesis/download"
# "/Users/manuel/development/thesis/staging"
# "/Users/manuel/development/thesis/overlap"

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


class EncodeOverlapper(Overlapper):
    def __init__(self, download_path, staging_path, overlap_path, storage_path):
        super(EncodeOverlapper, self).__init__(download_path, staging_path, overlap_path, storage_path)

    def overlap_bed_files(self, bed_filepath, bed_overlap_with_filepath, min_overlap):
        bed = BedTool(bed_filepath)
        intesect_bed = BedTool(bed_overlap_with_filepath)

        bed_intersect_bed = bed.intersect(intesect_bed, loj=True, f=min_overlap)
        return bed_intersect_bed

    def overlap_filtered_with_fantom_permissive(self, assembly='hg19', method=None, min_overlap=0.1, export_temp=False):
        storage_layer = HdfStoreManager(self.storage_path)

        encode_file_name = utils.build_filtered_file_name(assembly, method)
        encode_file_name_bed = utils.build_bed_file_name(encode_file_name)

        fantom_file_name = utils.build_permissive_file_name()
        fantom_file_name_bed = utils.build_bed_file_name(fantom_file_name)

        encode_bed = storage_layer.read_bed_file('encode_staging.hdf', encode_file_name_bed)
        fantom_bed = storage_layer.read_bed_file('fantom_staging.hdf', fantom_file_name_bed)

        print("ENCODE bed file from staging:", encode_file_name_bed)
        print("FANTOM bed file from staging:", fantom_file_name_bed)
        print("Starting overlap...")

        encode_intersect_fantom = encode_bed.intersect(fantom_bed, loj=True, f=min_overlap)
        print(len(encode_intersect_fantom), "ENCODE.intersect(FANTOM) left-outer-join results")

        overlap_column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'FA_chrom', 'FA_start', 'FA_end',
                                'FA_name', 'FA_score', 'FA_strand']

        overlap_df = encode_intersect_fantom.to_dataframe(names=overlap_column_names)

        overlap_df['size'] = overlap_df.apply(lambda row: self.compute_size(row), axis=1)

        if export_temp:
            output_filename_temp = encode_file_name + "_FANTOM_overlapped_temp.csv"
            print("Exporting temporary overlapped file (no ENCODE details merge) to temp:", output_filename_temp)
            storage_layer.store_dataframe_fixed(overlap_df, 'temp_overlaps.hdf', output_filename_temp)

        print("Adding details from FANTOM...")

        overlap_df['FA_size'] = overlap_df.apply(lambda row: self.compute_size(row, 'FA'), axis=1)
        overlap_df['FA_method'] = overlap_df.apply(lambda row: self.set_col_value(row, 'FA_name', 'CAGE_TCs'), axis=1)
        overlap_df['FA_ovlp_len'] = overlap_df.apply(lambda row: self.compute_ovlp_len(row, 'FA'), axis=1)
        overlap_df['FA_ovlp_pct'] = overlap_df.apply(lambda row: self.compute_ovlp_pct(row, 'FA'), axis=1)
        overlap_df['FA_encyclopedia'] = overlap_df.apply(lambda row: self.set_col_value(row, 'FA_name', 'FANTOM'),
                                                         axis=1)
        overlap_df['overlap_name'] = overlap_df.apply(lambda row: self.set_ovlp_name(row, 'name', 'FA_name'), axis=1)

        if export_temp:
            output_filename_temp = encode_file_name + "_FANTOM_overlapped_temp_detail.csv"
            print("Exporting temporary overlapped detailed file (no ENCODE details merge) to temp:",
                  output_filename_temp)
            storage_layer.store_dataframe_fixed(overlap_df, 'temp_overlaps.hdf', output_filename_temp)

        print("ENCODE details file:", encode_file_name)
        encode_details_df = storage_layer.read_dataframe('encode_staging.hdf', encode_file_name)

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
             'encyclopedia', 'overlap_name', 'FA_chrom', 'FA_start', 'FA_end',
             'FA_name', 'FA_score', 'FA_size', 'FA_method', 'FA_ovlp_len',
             'FA_ovlp_pct', 'FA_encyclopedia']]

        output_filename = encode_file_name + "_FANTOM_overlapped"

        print("Exporting overlaps file to overlaps:", output_filename,
              "- queryable columns are [biosample_term_name, biosample_type]")

        storage_layer.store_dataframe(overlap_full_detail_encode_df, 'encode_overlaps.hdf', output_filename,
                                      ['biosample_term_name', 'biosample_type'])

        output_bed_filename = output_filename + '_bed'
        print(
            "Exporting overlaps file in bed format "
            "(substituting names with overlaps name and considering only valid intersections) to overlaps:",
            output_bed_filename)

        overlap_full_detail_encode_df['name'] = overlap_full_detail_encode_df['overlap_name']
        full_bed_df_bed = overlap_full_detail_encode_df.query('FA_ovlp_pct > 0')[
            ['chrom', 'start', 'end', 'name', 'score', 'strand']]
        full_bed_df_bed.drop_duplicates(inplace=True)

        print("Storing", len(full_bed_df_bed), "unique overlaps...")
        storage_layer.store_dataframe_fixed(full_bed_df_bed, 'encode_overlaps.hdf', output_bed_filename)

        print("Completed")

    def overlap_filtered_with_dbsuper(self, assembly='hg19', method=None, min_overlap=0.1, export_temp=False):
        encode_file_name = utils.build_filtered_file_name(assembly, method)

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

    def overlap_filtered_with_roadmap(self, assembly='hg19', method=None, min_overlap=0.1, export_temp=False,
                                      use_test=False):
        encode_file_name = utils.build_filtered_file_name(assembly, method)

        encode_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".csv"
        roadmap_filepath = self.staging_path + "/EpigenomicsRoadmap/roadmap_metadata.csv"
        encode_bed_file_path = self.staging_path + "/ENCODE/filtered/" + encode_file_name + ".bed"
        roadmap_bed_filepath = self.staging_path + "/EpigenomicsRoadmap/processed.bed"

        if use_test:
            print("Using test bed files...")
            encode_bed_file_path = self.staging_path + "/ENCODE/filtered/filtered_test.bed"
            roadmap_bed_filepath = self.staging_path + "/EpigenomicsRoadmap/processed_test.bed"

        encode_bed = BedTool(encode_bed_file_path)
        roadmap_bed = BedTool(roadmap_bed_filepath)

        print("ENCODE bed file:", encode_bed_file_path)
        print("ROADMAP bed file:", roadmap_bed_filepath)
        print("Starting overlap...")

        encode_intersect_roadmap = encode_bed.intersect(roadmap_bed, loj=True, f=min_overlap)

        print(len(encode_intersect_roadmap), "ENCODE.intersect(ROADMAP) results")

        overlap_column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'RO_chrom', 'RO_start', 'RO_end',
                                'RO_name', 'RO_score', 'RO_strand']

        overlap_df = encode_intersect_roadmap.to_dataframe(names=overlap_column_names)
        overlap_df['size'] = overlap_df.apply(lambda row: self.compute_size(row), axis=1)

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_ROADMAP_overlapped_temp.csv"
            print("Exporting temporary overlapped file (no ENCODE details merge) to:", output_filename_temp)
            overlap_df.to_csv(output_filename_temp, index=None, sep='\t')

        print("ROADMAP details file:", roadmap_filepath)
        full_roadmap_details_df = pd.read_csv(roadmap_filepath, sep='\t')
        full_roadmap_details_df.columns = ['EID', 'RO_biosample_type', 'RO_biosample_group', 'RO_biosample_name',
                                           'RO_biosample_anatomy', 'imported']
        full_roadmap_details_df = full_roadmap_details_df.drop('imported', axis=1)

        # Adding ROADMAP EID to overlap_df
        overlap_df['EID'] = overlap_df.apply(lambda row: row['RO_name'][8:12], axis=1)

        if use_test:
            print(overlap_df)

        print("Merging details from ROADMAP...")
        overlap_full_detail_df = overlap_df.merge(
            full_roadmap_details_df[
                ['EID', 'RO_biosample_type', 'RO_biosample_group', 'RO_biosample_name', 'RO_biosample_anatomy']],
            how='left', on='EID')

        overlap_full_detail_df = overlap_full_detail_df.drop('EID', axis=1)

        overlap_full_detail_df['RO_size'] = overlap_full_detail_df.apply(lambda row: self.compute_size(row, 'RO'),
                                                                         axis=1)
        overlap_full_detail_df['RO_method'] = overlap_full_detail_df.apply(
            lambda row: self.set_col_value(row, 'RO_name', 'DNase'),
            axis=1)
        overlap_full_detail_df['RO_ovlp_len'] = overlap_full_detail_df.apply(
            lambda row: self.compute_ovlp_len(row, 'RO'), axis=1)
        overlap_full_detail_df['RO_ovlp_pct'] = overlap_full_detail_df.apply(
            lambda row: self.compute_ovlp_pct(row, 'RO'), axis=1)
        overlap_full_detail_df['RO_encyclopedia'] = overlap_full_detail_df.apply(
            lambda row: self.set_col_value(row, 'RO_name', 'ROADMAP'),
            axis=1)

        overlap_full_detail_df['RO_biosample_type'].fillna('.', inplace=True)
        overlap_full_detail_df['RO_biosample_group'].fillna('.', inplace=True)
        overlap_full_detail_df['RO_biosample_name'].fillna('.', inplace=True)
        overlap_full_detail_df['RO_biosample_anatomy'].fillna('.', inplace=True)

        if export_temp:
            output_filename_temp = self.overlap_path + "/" + encode_file_name + "_ROADMAP_overlapped_temp_detail.csv"
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
             'encyclopedia', 'RO_chrom', 'RO_start', 'RO_end',
             'RO_name', 'RO_score', 'RO_strand', 'RO_size', 'RO_method', 'RO_biosample_type', 'RO_biosample_group',
             'RO_biosample_name', 'RO_biosample_anatomy', 'RO_ovlp_len', 'RO_ovlp_pct', 'RO_encyclopedia']]

        output_filename = self.overlap_path + "/" + encode_file_name + "_ROADMAP_overlapped.csv"
        print("Exporting overlapped file to:", output_filename)
        overlap_full_detail_encode_df.to_csv(output_filename, index=None, sep='\t')

        print("Completed")
