import psutil
import time

from asyncio import sleep, subprocess
from collections import namedtuple
from pathlib import Path
from typing import List, Tuple

EvalResult = namedtuple(
    'EvalResult', ['tsv_path', 'runtime_s', 'peak_memory_MiB'])


async def _run_camel(
        work_dir: Path,
        proc_path: Path,
        arg_list: List[str]) -> Tuple[Path, float, float]:
    camel_dir = work_dir.joinpath('camel_out')
    camel_reads = camel_dir.joinpath('reads.fa')

    start, end = time.time(), time.time()
    peak_memory_usage = 0

    arg_list += ['-o', camel_dir.as_posix()]
    proc = await subprocess.create_subprocess_exec(proc_path, *arg_list)

    is_running = True
    while is_running:
        if proc.returncode:
            is_running = False
        else:
            try:
                pshandle = psutil.Process(proc.pid)

                end = time.time()
                peak_memory_usage = max(
                    peak_memory_usage,
                    pshandle.memory_info().rss)
            except:
                pass
        await sleep(0)

    return camel_reads, end - start, (peak_memory_usage / (2 ** 20))


async def _run_raven(
        work_dir: Path,
        reads_path: Path,
        threads: int) -> Path:
    raven_asm_path = work_dir.joinpath('raven_asm_path.fa')
    args = ['-t', str(threads), '-p', '0',
            '--disable-checkpoints', reads_path.as_posix()]

    with open(raven_asm_path, 'w') as asm_dst:
        proc = await subprocess.create_subprocess_exec(
            'raven', *args, stdout=asm_dst)
        await proc.wait()
    return raven_asm_path


async def _run_quast(
        work_dir: Path,
        asm_path: Path,
        ref_path: Path,
        threads: int) -> Path:
    quast_dir = work_dir.joinpath('quast_eval')
    quast_args = [
        '--fast', '-t', threads,
        '-r', ref_path.as_posix(),
        '-o', quast_dir.as_posix(),
        asm_path.as_posix()
    ]
    proc = await subprocess.create_subprocess_exec(
        'quast.py', *quast_args)
    await proc.wait()

    return quast_dir.joinpath('report.tsv')


async def eval_correction(
        work_dir: Path,
        ref_path: Path,
        threads: int,
        exe_path: Path,
        exe_args: List[str],
) -> EvalResult:
    reads, runtime_s, peak_memory_MiB = await _run_camel(work_dir, exe_path, exe_args)
    raven_asm_path = await _run_raven(work_dir, reads, threads)
    report_path = await _run_quast(work_dir, raven_asm_path, ref_path, threads)

    return EvalResult(report_path, runtime_s, peak_memory_MiB)
