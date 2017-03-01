import pandas as pd
from pybedtools import BedTool
import os


class HdfStoreManager:
    def __init__(self, storage_path):
        self.storage_path = storage_path

    def store_exists(self, store_name):
        store_name = self.storage_path + "/" + store_name
        return os.path.exists(store_name)

    def create_store(self, store_name):
        if not self.store_exists(store_name):
            store_name = self.storage_path + "/" + store_name
            store = pd.HDFStore(store_name)
            store.close()

    def table_exists(self, store_name, table):
        store_name = self.storage_path + "/" + store_name
        store = pd.HDFStore(store_name)
        exists = table in store
        store.close()
        return exists

    def store_data_from_csv(self, source, store_name, table_name, queryable_cols, sep='\t'):
        store_name = self.storage_path + "/" + store_name
        df = pd.read_csv(source, sep=sep)
        df.reset_index(level=0, inplace=True)

        df.to_hdf(store_name, table_name, mode='a', format='table', data_columns=queryable_cols)

    def store_dataframe(self, df, store_name, table_name, queryable_cols=None, append=False):
        store_name = self.storage_path + "/" + store_name
        df.to_hdf(store_name, table_name, mode='a', format='table', data_columns=queryable_cols, append=append)

    def store_dataframe_fixed(self, df, store_name, table_name, append=False):
        store_name = self.storage_path + "/" + store_name
        df.to_hdf(store_name, table_name, mode='a', format='fixed', append=append)

    def store_bed_file(self, bed, store_name, bed_name):
        store_name = self.storage_path + "/" + store_name
        bed.to_dataframe().to_hdf(store_name, bed_name, mode='a', format='fixed')

    def read_dataframe(self, store_name, table_name):
        store_name = self.storage_path + "/" + store_name
        return pd.read_hdf(store_name, table_name)

    def read_dataframe_by_query(self, store_name, table_name, query):
        store_name = self.storage_path + "/" + store_name
        return pd.read_hdf(store_name, table_name, where=query)

    def read_bed_file_as_dataframe(self, store_name, bed_name):
        store_name = self.storage_path + "/" + store_name
        store = pd.HDFStore(store_name)
        return store.select(bed_name, auto_close=True)

    def read_bed_file(self, store_name, bed_name):
        return BedTool().from_dataframe(self.read_bed_file_as_dataframe(store_name, bed_name))
