from overlaps_db.analysis.overlap_analyzer import OverlapAnalyzer
from overlaps_db.store.hdf_store_manager import HdfStoreManager

import overlaps_db.utils.data_process_utils as utils
from overlaps_db.utils.timer import Timer


class EncodeFantomAnalyzer(OverlapAnalyzer):
    def __init__(self, storage_path):
        super(EncodeFantomAnalyzer, self).__init__(storage_path)

    def perform_overlap_analysis_with_repeatmasker(self, assembly='hg19', method=None, overlap_intervals=10,
                                                   samples_num=20, biosample_type=None, repeat_family=None):
        pass

    def perform_reldist_analsys_with_repeatmasker(self, assembly='hg19', method=None,
                                                  biosample_type=None, repeat_family=None):
        pass
