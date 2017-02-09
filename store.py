#! /usr/bin/env python
# -*- coding: utf-8 -*-


import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.store.hdf_store import HdfStore

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    storage_dir = "/Users/manuel/development/thesis/storage"
    source_file = "/Users/manuel/development/thesis/overlap/filtered_hg19DNase_H3K27ac_dbSUPER_overlapped.csv"
    store = HdfStore(storage_path=storage_dir)
    store.store_data(source=source_file, store_name="encode_overlaps.hdf", table_name="encode_dbsuper",
                     queryable_cols=['biosample_term_name', 'SE_associated_gene'])


if __name__ == '__main__':
    cli.preamble("Store")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    received_args = parser.parse_args()

    sys.exit(main(received_args))
