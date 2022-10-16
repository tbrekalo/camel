import sqlite3

from typing import List, Union

_QUAST_TSV_TO_COLUMN_NAME = {
    '# contigs': 'n_contigs',
    'Largest contig': 'largest_contig',
    'N50': 'N50',
    'NG50': 'NG50',
    'NA50': 'NA50',
    'NGA50': 'NGA50',
    '# mismatches per 100 kbp': 'mismatches_per_100kbp',
    '# indels per 100 kbp': 'indels_per_100kbp',
    'Largest alignment': 'largest_alignment',
    'Total aligned length': 'total_aligned_length',
}

_QUAST_COLUMNS = [
    'n_contigs INTEGER',
    'largest_contig INTEGER',
    'N50 INTEGER',
    'NG50 INTEGER',
    'NA50 INTEGER',
    'NGA50 INTEGER',
    'mismatches_per_100kbp REAL',
    'indels_per_100kbp REAL',
    'largest_alignment INTEGER',
    'total_aligned_length INTEGER',
]

_RUNTIME_COLUMNS = [
    'runtime_s INTEGER',
    'peak_memory_MiB INTEGER',
    'timestamp DATETIME DEFAULT CURRENT_TIMESTAMP',
]

_TYPE_MAP = {
    int: 'INTEGER',
    float: 'REAL',
    str: 'TEXT',
}

StrongDict = dict[str, Union[int, float, str]]


def parse_quast_tsv(file):
    ret = {}
    with open(file) as tsv:
        for line in tsv.readlines():
            tsv_name, data = line.strip().rsplit('\t')
            db_name = _QUAST_TSV_TO_COLUMN_NAME.get(tsv_name)
            if db_name is not None:
                ret[db_name] = data

    return ret


def _flags_to_columns(exe_flags: StrongDict) -> List[str]:
    dst = []
    for k, v in exe_flags.items():
        dst.append(f'{k} {_TYPE_MAP[type(v)]}')

    return dst


def _create_table_sql(table_name: str, exe_flags: StrongDict) -> str:
    cols = ',\n'.join(_flags_to_columns(exe_flags) +
                      _RUNTIME_COLUMNS + _QUAST_COLUMNS)

    return f'CREATE TABLE IF NOT EXISTS {table_name}(\n{cols})'


def get_db_connection(db_path) -> sqlite3.Connection:
    ''' get cursor to sqlite3 db instance '''
    return sqlite3.connect(db_path)


def create_table_for_exe(
        db_connection: sqlite3.Connection,
        table_name: str,
        exe_flags: StrongDict) -> None:

    db_connection.cursor().execute(
        _create_table_sql(table_name, exe_flags))


def insert_dict_to_runs_table(
        db_connection: sqlite3.Connection,
        table_name: str,
        data: dict[str, str]):
    cursor = db_connection.cursor()

    columns_str = ', '.join(data.keys())
    values_str = ', '.join(map(lambda x: f'"{str(x)}"', data.values()))
    cmd_str = f'INSERT INTO {table_name} ({columns_str}) VALUES ({values_str})'

    cursor.execute(cmd_str)
    db_connection.commit()
