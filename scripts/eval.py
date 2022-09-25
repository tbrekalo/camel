import argparse
import psutil
import subprocess
import time

from pathlib import Path

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
                peek_mem_usage = max(peek_mem_usage, pshandle.memory_info().rss)
            except:
                pass
            time.sleep(.5)
    
    return camel_end - camel_start, peek_mem_usage

if __name__ == '__main__':
    parser = create_argparser()
    args = parser.parse_args()

    def extract_unary(name):
        ''' extract argument from parsed args '''
        return getattr(args, name)[0]

    work_dir = Path(extract_unary('work_dir'))
    camel_executable = extract_unary('executable')

    camel_args = [
        str(x) for k in ['threads', 'n_overlaps', 'window_length']
        for x in (f'--{k}', extract_unary(k))] + ['--dst', work_dir.as_posix()]
    camel_args.extend([
        getattr(args, 'overlaps'),
        getattr(args, 'reads'),
    ])

    elapsed, peek_mem = popen_time_mem(camel_args)
    print(f'peek mem usage : {peek_mem} MiB') 
