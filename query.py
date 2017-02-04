#! /usr/bin/env python
# -*- coding: utf-8 -*-


import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.query.parallel_load import ParallelLoad

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    file_name = "/Users/manuel/development/thesis/overlap/filtered_hg19DNase_H3K27ac_dbSUPER_overlapped.csv"
    parallel_load = ParallelLoad(file_name, args.processes, args.chunks_size)
    parallel_load.test()


if __name__ == '__main__':
    cli.preamble("Query")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument("--processes", action="store", dest="processes", type=int, default=1,
                        help='Number of processes, default is 1')
    parser.add_argument('--chunks_size', action="store", dest="chunks_size", type=int, default=100000,
                        help='Chunks size, default is 100,000')

    received_args = parser.parse_args()
    log.info("Processes: %d, Chunks size: %d", int(received_args.processes), int(received_args.chunks_size))

    sys.exit(main(received_args))
