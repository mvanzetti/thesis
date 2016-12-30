import json
import os
import urllib.request

import pandas as pd
import requests
from pandas.io.json import json_normalize


# Example:
# search_url = self.main_url + "/search/?type=Annotation&annotation_type=enhancer-like+regions&frame=object&limit=all"
# directory = "/Users/manuel/development/thesis/download/ENCODE"

class ENCODE:
    headers = {'accept': 'application/json'}
    main_url = "https://www.encodeproject.org"

    def __init__(self, type_param, annotation_type_param):
        self.search_url = \
            self.main_url + "/search/" + \
            "?type=" + type_param + \
            "&annotation_type=" + annotation_type_param + \
            "&frame=object&limit=all"

    def download(self, directory):

        response = requests.get(self.search_url, headers=self.headers)

        response_json_dict = response.json()

        graph = response_json_dict["@graph"]

        print("Saving to main dir ->", directory)
        if not os.path.exists(directory):
            print("Directory does not exist, creating...")
            os.makedirs(directory)

        file_list_name = directory + "/enhancer-like-annotations.csv"

        df = pd.DataFrame()

        if not os.path.exists(file_list_name):
            for i, element in enumerate(graph):
                temp_df = json_normalize(element)
                temp_df['index'] = i
                temp_df['imported'] = False
                df = df.append(temp_df)
                i += 1

            df.set_index('index')
            df.to_csv(file_list_name, sep='\t')

        else:
            df = pd.read_csv(file_list_name, sep='\t')

        count_downloaded = df.query('imported == True')['index'].count()
        print("Annotations entries to download:", count_downloaded, "of a total of", df['index'].count())

        for element in graph:
            print("Processing...", element["accession"], ":", element["description"])

            accession_directory = directory + "/" + element["accession"]

            if not os.path.exists(accession_directory):
                os.makedirs(accession_directory)

            file_path = accession_directory + "/info.txt"
            current_uuid = element["uuid"]
            element_index = df.query('uuid == @current_uuid').index

            if os.path.exists(file_path):
                print("info.txt exists, skipping")
            else:
                f = open(file_path, 'w')
                f.write(json.dumps(element, indent=4, separators=(',', ': ')))
                f.close()

            for file in element["files"]:
                info_file_url = self.main_url + file
                file_response = requests.get(info_file_url, headers=self.headers)
                file_response_dict = file_response.json()

                download_file_url = self.main_url + file_response_dict["href"]

                download_file_directory = accession_directory + file
                if not os.path.exists(download_file_directory):
                    os.makedirs(download_file_directory)

                filename = file_response_dict["href"].split("/")[-1]

                if "bed.gz" not in filename:
                    continue

                download_file_name = download_file_directory + filename

                if os.path.exists(download_file_name):
                    df = df.set_value(element_index, 'imported', True)
                    print("     File found, skipping ", filename)
                else:
                    print("     Downloading...", filename)
                    urllib.request.urlretrieve(download_file_url, filename=download_file_name)

            df = df.set_value(element_index, 'imported', True)
            df.to_csv(file_list_name, sep='\t')

            print("Completed")
