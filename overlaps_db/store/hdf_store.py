import pandas as pd


class HdfStore:
    def __init__(self, storage_path):
        self.storage_path = storage_path

    def store_data(self, source, store_name, table_name, queryable_cols):
        store_name = self.storage_path + "/" + store_name
        df = pd.read_csv(source, sep='\t')
        df.reset_index(level=0, inplace=True)

        df.to_hdf(store_name, table_name, mode='w', format='table', data_columns=queryable_cols)
