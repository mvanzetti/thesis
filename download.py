#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.data_process.roadmap_downloader import EpigenomicsRoadmapDownloader

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    download_dir = "/Users/manuel/development/thesis/download"
    staging_dir = "/Users/manuel/development/thesis/staging"
    overlap_dir = "/Users/manuel/development/thesis/overlap"

    if args.source == 'ROADMAP':
        log.info("Downloading Epigenomics Roadmap data")
        downloader = EpigenomicsRoadmapDownloader(download_directory=download_dir)
        downloader.download(args.force)

        log.info("Download completed")

    else:
        print("selected source not implemented yet")


if __name__ == '__main__':
    cli.preamble("Download")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument("source", help='Select one of the available sources',
                        choices=['ENCODE', 'ROADMAP', 'FANTOM', 'dbSUPER', 'All'])
    parser.add_argument('--force', action="store", dest="force",
                        help='If specified, force download even if resources are already present in the local dirs',
                        default=False)

    received_args = parser.parse_args()
    log.info("Downloading data for %s", received_args.source)
    log.info("Optional parameters: Force: %s", received_args.force)

    sys.exit(main(received_args))
