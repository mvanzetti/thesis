import urllib.request
import requests
import os
from bs4 import BeautifulSoup


class EpigenomicsRoadmapDownloader:
    headers = {'accept': 'application/json'}
    main_url = "http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/"

    def __init__(self):
        pass

    def download(self, directory):
        print("Sending request to", self.main_url, "...")
        response = requests.get(self.main_url, headers=self.headers)
        print("Response status is", response.status_code)

        html = response.text

        print("Parsing html response...")
        soup = BeautifulSoup(html, 'lxml')

        list_urls = soup.find_all('a')

        if not os.path.exists(directory):
            print("Directory does not exist, creating...")
            os.makedirs(directory)

        print("Saving to main dir ->", directory)
        for url in list_urls:
            file_name = url['href']
            newurl = self.main_url + file_name
            if newurl.endswith('bed.gz'):
                print("Downloading", file_name, "...")
                urllib.request.urlretrieve(newurl, directory + "/" + file_name)

        print("Completed")

# TODO Test only
roadmap_path = "/Users/manuel/development/thesis/download/EpigenomicsRoadmap"

downloader = EpigenomicsRoadmapDownloader()
downloader.download(directory=roadmap_path)
