import argparse
import logging
import sys

from overlaps_db.cli import defaults

log = logging.getLogger(sys.argv[0])


def connect_to_db(server=None):
    # imagine a connection to the file storage, etc
    pass


def start_logging(name=__name__, level=logging.DEBUG):
    logging.basicConfig(level=level, format="%(asctime)s : %(processName)s : %(levelname)s : %(name)s : %(message)s")
    return logging.getLogger(name)


def get_default_parser():
    parser = argparse.ArgumentParser()
    #    parser.add_argument('--user', action="store", dest="user", default='gpadmin')
    return parser


def print_banner(banner_path):
    try:
        f = open(banner_path, 'r')
        file_contents = f.read()
        print(file_contents)
        f.close()
    except FileNotFoundError as error:
        log.warning("Banner is missing:", error)
        print(defaults.banner_default_message)


def preamble(script_name=''):
    log.info("OverlapsDB " + script_name + " \U0001F680 \U0001F680 \U0001F680 \U0001F680 \U0001F680")
