import pandas as pd
import multiprocessing as mp
from overlaps_db.utils.timer import Timer


class ParallelLoad:
    def __init__(self, file_name, processes, chunk_size):
        self.file_name = file_name
        self.processes = processes
        self.chunk_size = chunk_size

    def process_frame(self, df):
        # process data frame
        return len(df)

    def test(self):
        timer = Timer()
        timer.start()
        reader = pd.read_table(self.file_name, chunksize=self.chunk_size)
        pool = mp.Pool(self.processes)

        funclist = []
        for df in reader:
            # process each data frame
            f = pool.apply_async(self.process_frame, [df])
            funclist.append(f)

        pool.close()
        print(timer.elapsed())
