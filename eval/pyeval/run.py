import asyncio
import psutil
import shutil
import subprocess
import time

from datetime import datetime
from pathlib import Path
from typing import Callable, Tuple

from db.orm.models.args import Args
from db.orm.models.benchmark import Benchmark
from db.orm.models.quast import Quast
from db.context import DBContext

from config import EvalConfig, StrongDict
from util import parse_quast_tsv


async def _run_correction(
        spawn: Callable[[], subprocess.Popen]) -> Tuple[float, float]:
    start, end = time.perf_counter(), time.perf_counter()
    peak_memory_usage = 0

    with spawn() as proc:
        pid = proc.pid
        pshandle = psutil.Process(pid)

        while proc.poll() is None:
            proc.returncode
            try:
                pshandle = psutil.Process(proc.pid)

                end = time.perf_counter()
                peak_memory_usage = max(
                    peak_memory_usage,
                    pshandle.memory_info().rss)
            except:
                pass
            await asyncio.sleep(1)

    return end - start, (peak_memory_usage / (2 ** 20))


async def _run_camel(
        work_dir: Path,
        executable: str,
        args: StrongDict,) -> Tuple[Path, float, float]:

    arg_list = []
    for k, v in args.items():
        arg_list.extend([f'--{k}', str(v)])

    camel_reads = work_dir.joinpath('camel_reads.fa')
    arg_list = arg_list + ['--out', work_dir.as_posix()]

    def spawn_fn():
        return subprocess.Popen([executable, *arg_list])

    runtime_s, peak_memory_MiB = await _run_correction(spawn_fn)
    return camel_reads, runtime_s, peak_memory_MiB


async def _run_racon(
        work_dir: Path,
        executable: str,
        args: StrongDict,) -> Tuple[Path, float, float]:
    racon_reads = work_dir.joinpath('racon_reads.fa')
    arg_list = [
        '-f',
        '--threads', str(args['threads']),
        '--window-length', str(args['window-length']),
        '--error-threshold', str(args['error-threshold']),
        args['reads'], args['overlaps'], args['reads']
    ]

    print(f'[pyeval::run_racon] args : {arg_list}')
    with open(racon_reads, 'w') as dst:
        def spawn_fn():
            return subprocess.Popen([executable, *arg_list],
                                    stdout=dst)

        runtime_s, peak_memory_MiB = await _run_correction(spawn_fn)
        return racon_reads, runtime_s, peak_memory_MiB


async def _run_raven(
        work_dir: Path,
        reads_path: Path,
        threads: int) -> Path:
    raven_asm_path = work_dir.joinpath('raven_asm_path.fa')
    args = ['-t', str(threads), '-p', '0',
            '--disable-checkpoints', reads_path.as_posix()]

    with open(raven_asm_path, 'w') as asm_dst:
        proc = await asyncio.subprocess.create_subprocess_exec(
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
        '--fast', '-t', str(threads),
        '-r', ref_path.as_posix(),
        '-o', quast_dir.as_posix(),
        asm_path.as_posix()
    ]
    proc = await asyncio.subprocess.create_subprocess_exec(
        'quast.py', *quast_args)
    await proc.wait()

    return quast_dir.joinpath('report.tsv')


def _get_version(
        executable: str) -> str:
    completion_handle = subprocess.run(
        [executable, '--version'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    return completion_handle.stderr.decode()


async def eval_correction(
        context: DBContext,
        cfg: EvalConfig,):
    timestamp = datetime.now().strftime('%d-%m-%Y_%H:%M')
    eval_dir = cfg.work_dir.joinpath(
        f'{cfg.executable.replace("/", "_")}_{timestamp}')
    if eval_dir.exists():
        shutil.rmtree(eval_dir)
    eval_dir.mkdir()

    try:
        context.executable = cfg.executable
        context.version = _get_version(cfg.executable)
        context.comment = cfg.comment

        context.args = Args(
            **{k.replace('-', '_'): v for k, v in cfg.args.items() if k != 'threads'})

        async def run_correction(eval_dir, executable, args):
            if executable == 'camel':
                return await _run_camel(eval_dir, 'camel', args)
            else:
                return await _run_racon(eval_dir, 'racon', args)

        reads, runtime_s, peak_memory_MiB = await run_correction(
            eval_dir,
            cfg.executable,
            cfg.args)

        context.benchmark = Benchmark(
            threads=cfg.args['threads'],
            runtime_s=runtime_s,
            peak_memory_mib=peak_memory_MiB)

        raven_asm_path = await _run_raven(eval_dir, reads, cfg.threads)
        report_path = await _run_quast(
            eval_dir, raven_asm_path, cfg.reference_path, cfg.threads)

        context.quast = Quast(**parse_quast_tsv(report_path))
    except Exception as err:
        print(err)
        raise err
    finally:
        shutil.rmtree(eval_dir)
