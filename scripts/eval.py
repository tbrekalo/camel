import argparse
import pprint
import psutil
import subprocess
import time

from pathlib import Path

QUAST_COLUMNS = [
    '# contigs',
    'Largest contig',
    'N50',
    'NG50',
    'NA50',
    'NGA50',
    '# mismatches per 100 kbp',
    '# indels per 100 kbp',
    'Largest alignment',
    'Total aligned length',
]


def parse_quast_tsv(file):
    ret = {}
    with open(file) as tsv:
        for line in tsv.readlines():
            name, data = line.strip().rsplit('\t')
            if name in QUAST_COLUMNS:
                ret[name] = data

    return ret


def create_argparser():
    ''' create argparser with arguments '''
    dst = argparse.ArgumentParser()

    def add_unary_arg(*args, **kwargs):
        dst.add_argument(*args, nargs=1, action='store', **kwargs)

    add_unary_arg('-c', '--comment', type=str,  required=True)
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


def popen_time_mem(arg_list):
    camel_start, camel_end = time.time(), time.time()
    peek_mem_usage = 0
    with subprocess.Popen([camel_executable, *camel_args]) as camel_proc:
        pid = camel_proc.pid
        pshandle = psutil.Process(pid)
        while camel_proc.poll() is None:
            try:
                camel_end = time.time()
                peek_mem_usage = max(
                    peek_mem_usage, pshandle.memory_info().rss)
            except:
                pass
            time.sleep(.5)

    return camel_end - camel_start, peek_mem_usage


if __name__ == '__main__':
    parser = create_argparser()
    args = parser.parse_args()
    log_data = {}

    def extract_unary_as_str(name):
        ''' extract argument from parsed args '''
        return str(getattr(args, name)[0])

    work_dir = Path(extract_unary_as_str('work_dir'))
    camel_executable = extract_unary_as_str('executable')

    camel_args = [
        x for k in ['threads', 'n_overlaps', 'window_length']
        for x in (f'--{k}', extract_unary_as_str(k))] + ['--dst', work_dir.as_posix()]
    camel_args.extend([
        getattr(args, 'overlaps'),
        getattr(args, 'reads'),
    ])

    runtime, peak_memory = popen_time_mem(camel_args)
    log_data = {
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

        pprint.pprint(log_data, indent=2, width=2)
