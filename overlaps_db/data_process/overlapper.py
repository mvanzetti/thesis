class Overlapper:
    def __init__(self, download_path, staging_path, overlap_path):
        self.download_path = download_path
        self.staging_path = staging_path
        self.overlap_path = overlap_path

    @staticmethod
    def set_col_value(row, col_name, value):
        if row[col_name] == '.':
            return '.'
        return value

    @staticmethod
    def compute_size(row, prefix=None):
        col_name = prefix + '_name' if prefix else 'name'
        if row[col_name] == '.':
            return 0
        col_end = prefix + '_end' if prefix else 'end'
        col_start = prefix + '_start' if prefix else 'start'
        size = abs(row[col_end] - row[col_start])
        return size

    @staticmethod
    def compute_ovlp_len(row, prefix):
        col_name = prefix + '_name'
        col_end = prefix + '_end'
        col_start = prefix + '_start'
        if row[col_name] == '.':
            return 0
        value = min(row['end'], row[col_end]) - max(row['start'], row[col_start])
        return value

    @staticmethod
    def compute_ovlp_pct(row, prefix):
        col_len = prefix + '_ovlp_len'
        size = abs(row['end'] - row['start'])
        value = row[col_len] / size * 100.0
        return value
