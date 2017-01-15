import os

import pandas as pd
import urllib.request
import requests
import shutil
from bs4 import BeautifulSoup


class DbSuperDownloader:
    headers = {'accept': 'application/json'}
    main_url = "http://bioinfo.au.tsinghua.edu.cn/dbsuper/"

    def __init__(self, genome_type, cell_type='all'):
        self.search_url = self.main_url + "adv_search.php?genome=" + genome_type
        self.search_bed_file_url = self.main_url + "data/bed/" + genome_type + "/all_" + genome_type + "_bed.bed"

        if cell_type == 'all':
            print("Selected cell type", cell_type)
            self.search_url += "&cell_type%5B%5D=" + cell_type + "&gene=&locus=&method=&submit=1"
        else:
            print("Selected all cell types")
            self.search_url += "&cell_type%5B%5D=C_001&cell_type%5B%5D=C_002&cell_type%5B%5D=C_003\
        &cell_type%5B%5D=C_004&cell_type%5B%5D=C_005&cell_type%5B%5D=C_006&cell_type%5B%5D=C_007\
        &cell_type%5B%5D=C_008&cell_type%5B%5D=C_009&cell_type%5B%5D=C_010&cell_type%5B%5D=C_011\
        &cell_type%5B%5D=C_012&cell_type%5B%5D=C_013&cell_type%5B%5D=C_014&cell_type%5B%5D=C_015\
        &cell_type%5B%5D=C_016&cell_type%5B%5D=C_017&cell_type%5B%5D=C_018&cell_type%5B%5D=C_019\
        &cell_type%5B%5D=C_020&cell_type%5B%5D=C_021&cell_type%5B%5D=C_022&cell_type%5B%5D=C_023\
        &cell_type%5B%5D=C_024&cell_type%5B%5D=C_025&cell_type%5B%5D=C_026&cell_type%5B%5D=C_027\
        &cell_type%5B%5D=C_028&cell_type%5B%5D=C_029&cell_type%5B%5D=C_030&cell_type%5B%5D=C_031\
        &cell_type%5B%5D=C_032&cell_type%5B%5D=C_033&cell_type%5B%5D=C_034&cell_type%5B%5D=C_035\
        &cell_type%5B%5D=C_036&cell_type%5B%5D=C_037&cell_type%5B%5D=C_038&cell_type%5B%5D=C_039\
        &cell_type%5B%5D=C_096&cell_type%5B%5D=C_040&cell_type%5B%5D=C_041&cell_type%5B%5D=C_042\
        &cell_type%5B%5D=C_043&cell_type%5B%5D=C_044&cell_type%5B%5D=C_045&cell_type%5B%5D=C_046\
        &cell_type%5B%5D=C_047&cell_type%5B%5D=C_100&cell_type%5B%5D=C_048&cell_type%5B%5D=C_049\
        &cell_type%5B%5D=C_050&cell_type%5B%5D=C_128&cell_type%5B%5D=C_097&cell_type%5B%5D=C_051\
        &cell_type%5B%5D=C_052&cell_type%5B%5D=C_053&cell_type%5B%5D=C_054&cell_type%5B%5D=C_055\
        &cell_type%5B%5D=C_103&cell_type%5B%5D=C_056&cell_type%5B%5D=C_057&cell_type%5B%5D=C_058\
        &cell_type%5B%5D=C_059&cell_type%5B%5D=C_060&cell_type%5B%5D=C_061&cell_type%5B%5D=C_062\
        &cell_type%5B%5D=C_063&cell_type%5B%5D=C_093&cell_type%5B%5D=C_094&cell_type%5B%5D=C_095\
        &cell_type%5B%5D=C_064&cell_type%5B%5D=C_065&cell_type%5B%5D=C_102&cell_type%5B%5D=C_101\
        &cell_type%5B%5D=C_066&cell_type%5B%5D=C_104&cell_type%5B%5D=C_067&cell_type%5B%5D=C_068\
        &cell_type%5B%5D=C_069&cell_type%5B%5D=C_070&cell_type%5B%5D=C_071&cell_type%5B%5D=C_105\
        &cell_type%5B%5D=C_072&cell_type%5B%5D=C_073&cell_type%5B%5D=C_074&cell_type%5B%5D=C_075\
        &cell_type%5B%5D=C_076&cell_type%5B%5D=C_077&cell_type%5B%5D=C_078&cell_type%5B%5D=C_079\
        &cell_type%5B%5D=C_080&cell_type%5B%5D=C_081&cell_type%5B%5D=C_108&cell_type%5B%5D=C_109\
        &cell_type%5B%5D=C_082&cell_type%5B%5D=C_098&cell_type%5B%5D=C_099&cell_type%5B%5D=C_083\
        &cell_type%5B%5D=C_084&cell_type%5B%5D=C_085&cell_type%5B%5D=C_086&gene=&locus=&method=&submit=1"

    def download_annotations(self, directory, annotation_filename):
        print("Sending request to", self.main_url, "...")
        response = requests.get(self.search_url, headers=self.headers)
        print("Response status is", response.status_code)

        html = response.text

        print("Parsing html response...")
        soup = BeautifulSoup(html, 'lxml')
        table = soup.find_all('table')[0]

        htr = table.find_all('tr')[0]
        htds = htr.find_all('th')[0:10]
        columns = [elem.text for elem in htds]

        rows = []
        for tr in table.find_all('tr')[1:]:
            tds = tr.find_all('td')[0:10]
            rows.append([elem.text.strip() for elem in tds])

        df = pd.DataFrame.from_records(rows, columns=columns)
        df.index.name = 'index'

        print("Parsed", len(df), " rows with columns:", columns)

        print("Saving to main dir ->", directory)
        if not os.path.exists(directory):
            print("Directory does not exist, creating...")
            os.makedirs(directory)

        file_list_name = directory + "/" + annotation_filename

        df.to_csv(file_list_name, sep='\t')

        print("Annotation file saved to", file_list_name)

        print("Completed")

    def download_bed(self, directory, bed_filename=""):
        print("Downloading the full bed file from", self.search_bed_file_url)

        if not bed_filename:
            bed_file_name = self.search_bed_file_url.split("/")[-1]

        print("Saving to main dir ->", directory)
        if not os.path.exists(directory):
            print("Directory does not exist, creating...")
            os.makedirs(directory)

        bed_file_path = directory + "/" + bed_file_name
        with urllib.request.urlopen(self.search_bed_file_url) as response, open(bed_file_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print("Bed file saved to", bed_file_path)
        print("Completed")


# TODO Test only
dbsuper_path = "/Users/manuel/development/thesis/download/dbSUPER"
ann_filename = "super-enhancers-annotations.csv"

downloader = DbSuperDownloader(genome_type="hg19")
# downloader.download_annotations(directory=dbsuper_path, annotation_filename=ann_filename)
downloader.download_bed(directory=dbsuper_path)
