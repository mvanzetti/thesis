#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.data_process.encode_overlapper import EncodeOverlapper

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    download_dir = "/Users/manuel/development/thesis/download"
    staging_dir = "/Users/manuel/development/thesis/staging"
    overlap_dir = "/Users/manuel/development/thesis/overlap"

    if args.source == 'ENCODE':
        log.info("Initializing EncodeOverlapper")
        overlapper = EncodeOverlapper(download_dir, staging_dir, overlap_dir)

        if args.encyclopedia == 'FANTOM':
            log.info("Overlapping with FANTOM")
            overlapper.overlap_filtered_with_fantom_permissive(assembly=args.assembly, method=args.method,
                                                               min_overlap=args.min_overlap,
                                                               export_temp=args.export_temp)
        elif args.encyclopedia == 'dbSUPER':
            log.info("Overlapping with dbSUPER")
            overlapper.overlap_filtered_with_dbsuper(assembly=args.assembly, method=args.method,
                                                     min_overlap=args.min_overlap,
                                                     export_temp=args.export_temp)

        elif args.encyclopedia == 'ROADMAP':
            log.info("Overlapping with Epigenomics Roadmap")
            overlapper.overlap_filtered_with_roadmap(assembly=args.assembly, method=args.method,
                                                     min_overlap=args.min_overlap,
                                                     export_temp=args.export_temp,
                                                     use_test=args.use_test)

        log.info("Overlapping completed")
    else:
        print("selected source not implemented yet")


if __name__ == '__main__':
    cli.preamble("Analyze")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument("source", help='Select one of the available sources', choices=['ENCODE', 'ENCODE_FANTOM'])
    parser.add_argument("encyclopedia",
                        help='Select the encyclopedia to consider',
                        choices=['FANTOM', 'dbSUPER', 'ROADMAP'])
    parser.add_argument('--assembly', action="store", dest="assembly",
                        help='The assembly to use to build the overlaps files. \
                        For the human genome assembly, type hg19.',
                        default='hg19')
    parser.add_argument('--method', action="store", dest="method",
                        help='Filter the source enhancer list by the provided method',
                        default='DNase_H3K27ac')
    parser.add_argument('--minoverlap', action="store", dest="min_overlap",
                        help='The minimum overlap requested while overlapping.', default=0.1)
    parser.add_argument('--temp', action="store", dest="export_temp",
                        help='If specified, export temporary files representing different overlapping phases',
                        default=False)
    parser.add_argument('--test', action="store", dest="use_test",
                        help='If specified, consider test file while overlapping. Use only for test purposes',
                        default=False)


    received_args = parser.parse_args()
    log.info("Building OverlapsDB for %s in %s", received_args.source, received_args.encyclopedia)
    log.info("Optional parameters: Assembly: %s, Method: %s, Min Overlap: %s, Export temp: %s", received_args.assembly,
             received_args.method, received_args.min_overlap, received_args.export_temp)

    sys.exit(main(received_args))
