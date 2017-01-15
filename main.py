#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.data_process.encode_overlapper import EncodeOverlapper

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    operation = args.operation

    if operation == 'download':
        pass
    elif operation == 'process':
        pass
    elif operation == 'overlap':
        overlap(args)
    elif operation == 'hello':
        print("Hello!")


def init(args):
    # todo init dirs
    pass


def overlap(args):
    download_dir = "/Users/manuel/development/thesis/download"
    staging_dir = "/Users/manuel/development/thesis/staging"
    overlap_dir = "/Users/manuel/development/thesis/overlap"

    overlapper = EncodeOverlapper(download_dir, staging_dir, overlap_dir)
    overlapper.overlap_filtered_with_dbsuper(assembly=args.assembly, method=args.method, min_overlap=args.min_overlap,
                                             export_temp=args.export_temp)


if __name__ == '__main__':
    cli.preamble("Main")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument('-op', action="store", dest="operation",
                        help='Available operations: download, process, overlap', default='hello')
    parser.add_argument('-a', action="store", dest="assembly", help='Assembly', default='hg19')
    parser.add_argument('-m', action="store", dest="method", help='Method', default='DNase_H3K27ac')
    parser.add_argument('-o', action="store", dest="min_overlap", help='Minimum Overlap', default=0.1)
    parser.add_argument('-e', action="store", dest="export_temp", help='Export Temporary Overlaps', default=False)

    received_args = parser.parse_args()
    sys.exit(main(received_args))
