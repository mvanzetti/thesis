#!/usr/bin/python
from data_process import encode_downloader
from data_process import encode_processor

encode_path = "/Users/manuel/development/thesis/download/ENCODE"
annotation_filename = "enhancer-like-annotations.csv"

print("Checking for ENCODE resources...")
downloader = encode_downloader.EncodeDownloader(type_param="Annotation", annotation_type_param="enhancer-like+regions")
downloader.download(directory=encode_path, annotation_filename=annotation_filename)

# print("Merging ENCODE resources...")
# output_filename = "enhancer-like-bed_refs.csv"
#
# merger = merger.EncodeBedMerger(encode_path)
# merger.prepare(annotation_filename, output_filename)

#blah