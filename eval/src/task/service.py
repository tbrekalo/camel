from pathlib import Path
from time import perf_counter
from typing import List

from psutil import Popen

from task.schemas import TaskConfig, TaskRun


def format_args(task_cfg: TaskConfig) -> List[str]:
    return [
        '--threads', str(task_cfg.threads),
        '--window-length', str(task_cfg.window_length),
        '--error-threshold', str(task_cfg.error_threshold),
    ]


def format_camel_args(task_cfg: TaskConfig) -> List[str]:
    return [
        task_cfg.exe,
        *format_args(task_cfg),
        str(task_cfg.reads_path),
        str(task_cfg.overlaps_path),
    ]


def format_racon_args(task_cfg: TaskConfig) -> List[str]:
    return format_camel_args(task_cfg) + [
        str(task_cfg.reads_path),
        '-f'
    ]


def create_spawn_list(task_cfg: TaskConfig) -> List[str]:
    if task_cfg.exe.endswith('camel'):
        return format_camel_args(task_cfg)
    return format_racon_args(task_cfg)


def monitored_run(output_dir: Path, task_cfg: TaskConfig) -> TaskRun:
    reads_path = output_dir.joinpath('reads.fa')
    with open(reads_path, 'w+') as f:
        with Popen(create_spawn_list(task_cfg), stdout=f) as proc:
            peak_memory = 0
            time_begin = time_end = perf_counter()

            while proc.poll() is None:
                curr_mem = proc.memory_info().rss
                time_end = perf_counter()

                if curr_mem is not None and curr_mem > peak_memory:
                    peak_memory = curr_mem

    ret = TaskRun(
        peak_memory_mib=peak_memory / (2 ** 20),
        runtime_s=time_end - time_begin,
    )

    return ret
