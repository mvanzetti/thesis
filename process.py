#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.data_process.roadmap_processor import EpigenomicsRoadmapProcessor
from overlaps_db.data_process.encode_processor import EncodeBedProcessor
from overlaps_db.data_process.fantom_processor import FantomBedProcessor
from overlaps_db.data_process.repeatmasker_processor import RepeatMaskerProcessor

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    download_dir = "/Users/manuel/development/thesis/download"
    staging_dir = "/Users/manuel/development/thesis/staging"
    storage_dir = "/Users/manuel/development/thesis/storage"
    overlap_dir = "/Users/manuel/development/thesis/overlap"

    if args.source == 'ENCODE':
        log.info("Processing ENCODE data")
        processor = EncodeBedProcessor(download_path=download_dir, staging_path=staging_dir, storage_path=storage_dir)
        processor.prepare()
        processor.filter(assembly=args.assembly, method=args.method, force=args.force)
        log.info("Process completed")

    elif args.source == 'FANTOM':
        log.info("Processing FANTOM data")
        processor = FantomBedProcessor(download_path=download_dir, staging_path=staging_dir, storage_path=storage_dir)
        processor.process_permissive()
        log.info("Process completed")

    elif args.source == 'ROADMAP':
        log.info("Processing Epigenomics Roadmap data")
        processor = EpigenomicsRoadmapProcessor(download_path=download_dir, staging_path=staging_dir,
                                                storage_path=storage_dir)
        processor.process()
        log.info("Process completed")

    elif args.source == 'RepeatMasker':
        log.info("Processing Repeat Masker data")
        processor = RepeatMaskerProcessor(download_path=download_dir, staging_path=staging_dir,
                                          storage_path=storage_dir)
        processor.process(assembly=args.assembly)
        log.info("Process completed")

    else:
        print("selected source not implemented yet")


if __name__ == '__main__':
    cli.preamble("Process")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument("source", help='Select one of the available sources',
                        choices=['ENCODE', 'ROADMAP', 'FANTOM', 'dbSUPER', 'RepeatMasker'])
    parser.add_argument('--assembly', action="store", dest="assembly",
                        help='The assembly to use to filter and process files. \
                        For the human genome assembly, type hg19.',
                        default='hg19')
    parser.add_argument('--method', action="store", dest="method",
                        help='When available, filter the source list by the provided method',
                        default='DNase_H3K27ac')
    parser.add_argument('--force', action="store", dest="force",
                        help='If specified, force process even if resources are already present in the local dirs',
                        default=False)

    received_args = parser.parse_args()
    log.info("Processing data for %s", received_args.source)
    log.info("Optional parameters: Assembly: %s, Method: %s, Force: %s", received_args.assembly, received_args.method,
             received_args.force)

    sys.exit(main(received_args))
