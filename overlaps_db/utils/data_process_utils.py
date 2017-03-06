def join_col_values(row, col, sep):
    return "[" + sep.join(row[col]) + "]"


def stringify_columns(df, columns, sep=','):
    for index, row in df.iterrows():
        for col in columns:
            df.set_value(index, col, join_col_values(row, col, sep))


def build_filtered_file_name(assembly, method):
    encode_file_name = "filtered"
    if assembly:
        encode_file_name += "_" + assembly
    if method:
        encode_file_name += "_" + method
    return encode_file_name


def build_biosample_file_name(biosample, assembly, method):
    file_name = biosample.replace(" ", "_").replace(":", "_")
    if assembly:
        file_name += "_" + assembly
    if method:
        file_name += "_" + method
    return file_name


def build_permissive_file_name():
    return "permissive"


def build_repeatmasker_file_name(repeat_class):
    return repeat_class.replace('/', '_').replace('?', '_qm').replace('-', '_')


def build_bed_file_name(file_name):
    return file_name + "_bed"


def build_biosample_bed_file_name(biosample, assembly, method):
    return build_bed_file_name(build_biosample_file_name(biosample, assembly, method))
