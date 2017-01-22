import urllib.request
import requests
import os
import pandas as pd
from bs4 import BeautifulSoup


class EpigenomicsRoadmapDownloader:
    headers = {'accept': 'application/json'}
    main_url = "http://egg2.wustl.edu/roadmap/data/"

    def __init__(self, download_directory):
        self.download_directory = download_directory + "/EpigenomicsRoadmap"

    @staticmethod
    def make_directory(directory):
        if not os.path.exists(directory):
            print("Directory", directory, "does not exist, creating...")
            os.makedirs(directory)

    def download_metadata(self, metadata_filepath):
        self.make_directory(self.download_directory)

        download_url = self.main_url + "byFileType/metadata/EID_metadata.tab"
        file_path = self.download_directory + "/" + "EID_metadata.tab"

        urllib.request.urlretrieve(download_url, file_path)

        metadata_df = pd.read_csv(file_path, sep='\t')
        metadata_df = metadata_df[['EID', 'TYPE', 'GROUP', 'STD_NAME', 'ANATOMY']]
        metadata_df.columns = ['EID', 'biosample_type', 'biosample_group', 'biosample_name', 'biosample_anatomy']
        metadata_df['imported'] = False

        metadata_df.to_csv(metadata_filepath, sep='\t', index=False)

    def download(self, force=False):
        metadata_filename = self.download_directory + "/" + "roadmap_metadata.csv"

        if not os.path.exists(metadata_filename) and not force:
            self.download_metadata(metadata_filename)

        metadata_df = pd.read_csv(metadata_filename, sep='\t')

        to_download = len(metadata_df.query('imported == False'))
        print("Annotations entries to download:", to_download, "of a total of", len(metadata_df))

        download_url = self.main_url + "byDataType/dnase/BED_files_enh/"
        print("Sending request to", download_url, "...")
        response = requests.get(download_url, headers=self.headers)
        print("Response status is", response.status_code)

        html = response.text

        print("Parsing html response...")
        soup = BeautifulSoup(html, 'lxml')

        list_urls = soup.find_all('a')

        self.make_directory(self.download_directory)

        print("Saving to main dir ->", self.download_directory)
        for url in list_urls:
            file_name = url['href']
            eid = file_name[12:16]
            eid_df = metadata_df.query('EID == @eid')
            if len(eid_df) == 1:
                eid_index = eid_df.index[0]
                if not metadata_df.get_value(eid_index, 'imported'):
                    newurl = download_url + file_name
                    if newurl.endswith('bed.gz'):
                        print("Downloading", file_name, "...")
                        urllib.request.urlretrieve(newurl, self.download_directory + "/" + file_name)
                        metadata_df = metadata_df.set_value(eid_index, 'imported', True)
            else:
                print("Too metadata rows for", file_name)

        metadata_df.to_csv(metadata_filename, sep='\t', index=False)
        print("Completed")
