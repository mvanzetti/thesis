def join_col_values(row, col, sep):
    return "[" + sep.join(row[col]) + "]"


def stringify_columns(df, columns, sep=','):
    for index, row in df.iterrows():
        for col in columns:
            df.set_value(index, col, join_col_values(row, col, sep))
