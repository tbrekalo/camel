import argparse
import psutil
import subprocess
import time
import sqlite3

from pathlib import Path

QUAST_COLUMNS = {
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

CREATE_CAMEL_RUNS = '''
CREATE TABLE IF NOT EXISTS camel_runs(
    camel_version TEXT,
    comment TEXT,
    window_length INTEGER,
    n_overlaps INTEGER,
    threads INTEGER,
    reads_path TEXT,
    overlaps_path TEXT,
    reference_path TEXT,
    runtime_s INTEGER,
    peak_memory_MiB INTEGER,
    n_contigs INTEGER,
    largest_contig INTEGER,
    N50 INTEGER,
    NG50 INTEGER,
    NA50 INTEGER,
    NGA50 INTEGER,
    mismatches_per_100kbp REAL,
    indels_per_100kbp REAL,
    largest_alignment INTEGER,
    total_aligned_length INTEGER,
    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
);'''


def get_db_connection(db_path):
    ''' get cursor to sqlite3 db instance '''
    return sqlite3.connect(db_path)


def create_runs_table_if_not_exists(db_connection):
    return db_connection.cursor().execute(CREATE_CAMEL_RUNS)


def insert_dict_to_runs_table(db_connection, data):
    cursor = db_connection.cursor()

    columns_str = ', '.join(data.keys())
    values_str = ', '.join(map(lambda x: f'"{str(x)}"', data.values()))
    cmd_str = f'INSERT INTO camel_runs ({columns_str}) VALUES ({values_str})'

    cursor.execute(cmd_str)
    db_connection.commit()


def parse_quast_tsv(file):
    ret = {}
    with open(file) as tsv:
        for line in tsv.readlines():
            tsv_name, data = line.strip().rsplit('\t')
            db_name = QUAST_COLUMNS.get(tsv_name)
            if db_name is not None:
                ret[db_name] = data

    return ret


def create_argparser():
    ''' create argparser with arguments '''
    dst = argparse.ArgumentParser()

    def add_unary_arg(*args, **kwargs):
        dst.add_argument(*args, nargs=1, action='store', **kwargs)

    add_unary_arg('-d', '--database', type=str, default=['camel_dev.db'])
    add_unary_arg('-c', '--comment', type=str,  required=True)
    add_unary_arg('-v', '--version', type=str, required=True)
    add_unary_arg('-o', '--n_overlaps', type=int,  default=[480])
    add_unary_arg('-t', '--threads', type=int,  default=[1])
    add_unary_arg('-e', '--executable', type=str,
                  default=['./build/bin/camel_exe'])
    add_unary_arg('-r', '--reference', type=str, required=True)
    add_unary_arg('-w', '--window_length', type=int,  default=[500])
    add_unary_arg('--work_dir', type=str, default=['/tmp/camel_eval'])

    dst.add_argument('overlaps')
    dst.add_argument('reads')

    return dst


def popen_time_mem(proc_path, arg_list):
    start, end = time.time(), time.time()
    peek_mem_usage = 0
    with subprocess.Popen([proc_path, *arg_list]) as proc:
        pid = proc.pid
        pshandle = psutil.Process(pid)
        while proc.poll() is None:
            try:
                end = time.time()
                peek_mem_usage = max(
                    peek_mem_usage, pshandle.memory_info().rss)
            except:
                pass
            time.sleep(.5)

    return end - start, peek_mem_usage


if __name__ == '__main__':
    parser = create_argparser()
    args = parser.parse_args()
    log_data = {}

    def extract_unary_as_str(name):
        ''' extract argument from parsed args '''
        return str(getattr(args, name)[0])

    db_path = extract_unary_as_str('database')
    work_dir = Path(extract_unary_as_str('work_dir'))

    camel_version = extract_unary_as_str('version')
    camel_executable = extract_unary_as_str('executable')

    comment = extract_unary_as_str('comment')

    camel_args = [
        x for k in ['threads', 'n_overlaps', 'window_length']
        for x in (f'--{k}', extract_unary_as_str(k))] + ['--dst', work_dir.as_posix()]
    camel_args.extend([
        getattr(args, 'overlaps'),
        getattr(args, 'reads'),
    ])

    db_connection = get_db_connection(db_path)
    create_runs_table_if_not_exists(db_connection)
    runtime, peak_memory = popen_time_mem(camel_executable, camel_args)
    log_data = {
        'camel_version': camel_version,
        'comment': comment,
        'window_length': extract_unary_as_str('window_length'),
        'n_overlaps': extract_unary_as_str('n_overlaps'),
        'threads': extract_unary_as_str('threads'),
        'reads_path': getattr(args, 'reads'),
        'overlaps_path': getattr(args, 'overlaps'),
        'reference_path': extract_unary_as_str('reference'),
        'runtime_s': runtime,
        'peak_memory_MiB': (peak_memory / (2 ** 20)),
    }

    corrected_reads_path = work_dir.joinpath('reads.fa')
    raven_args = ['-t', extract_unary_as_str('threads'), '-p', '0', '--disable-checkpoints',
                  corrected_reads_path.as_posix()]

    raven_asm_path = work_dir.joinpath('raven_camel_asm.fa')
    with open(raven_asm_path, 'w') as raven_asm_dst:
        subprocess.run(
            ' '.join(['raven', *raven_args]), shell=True, stdout=raven_asm_dst)

        quast_results_dir = work_dir.joinpath('camel_quast')
        quast_args = [
            '-t', extract_unary_as_str('threads'), '--fast',
            '-r', extract_unary_as_str(
                'reference'), raven_asm_path.as_posix(),
            '-o', quast_results_dir.as_posix()]
        subprocess.run(' '.join(['quast.py', *quast_args]),
                       shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        quast_tsv_report = quast_results_dir.joinpath('report.tsv')
        log_data |= parse_quast_tsv(quast_tsv_report)

        insert_dict_to_runs_table(db_connection, log_data)
