class EncodeOverlapper:
    def __init__(self, staging_path, overlap_path):
        self.staging_path = staging_path
        self.overlap_path = overlap_path

    def overlapFilteredWithdbSuper(self, assembly=None, method=None, force=False):
        encode_file_name = "filtered_"
        if assembly:
            encode_file_name += assembly
        if method:
            encode_file_name += method

        encode_file_path = self.staging_path + "/filtered/" + encode_file_name + ".csv"

        # create a big bed file from the dataframe

