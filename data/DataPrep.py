#!/usr/bin/python
# from .EncodeDownloader import EncodeDownloader
from data import ResourceDownloader

print("Checking for ENCODE resources...")
downloader = ResourceDownloader.ENCODE("Annotation", "enhancer-like+regions")
downloader.download("/Users/manuel/development/thesis/download/ENCODE")

print("Merging ENCODE resources...")
# TODO

#blah