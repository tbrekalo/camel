import argparse
import json
import pathlib
import sys
from datetime import datetime

from exceptions import ReferenceNotFound
from task.schemas import TaskConfig, TaskInfo
from task.service import monitored_run

parser = argparse.ArgumentParser(
    prog='evaluation',
    description='evaluate read correction tool',
)

parser.add_argument(
    'config',
    help='utf-8 json file with tool configuration',
    type=str,
)

parser.add_argument(
    '-r', '--reference',
    help='reference reads',
    type=str,
    required=True,
)

parser.add_argument(
    '-o', '--output',
    help='output folder containing runtime information',
    type=str,
    required=True,
)

try:
    args = parser.parse_args()
    task_cfg = TaskConfig.parse_file(args.config)

    ref_path = pathlib.Path(args.reference)
    if not ref_path.exists():
        raise ReferenceNotFound(
            args.reference,
            "could not find reference on the disk"
        )

    output_dir = pathlib.Path(args.output)
    if not output_dir.exists():
        print('creating output dir', file=sys.stderr)
        output_dir.mkdir()

    runtime_dir_path = output_dir.joinpath(
        '{}_{}'.format(
            task_cfg.exe,
            datetime.now().strftime('%y-%m-%d_%H-%M')
        )
    )

    runtime_dir_path.mkdir()
    with open(runtime_dir_path.joinpath('info.json'), 'w+') as f:
        f.write(TaskInfo(
            task_config=task_cfg,
            task_run=monitored_run(runtime_dir_path, task_cfg),
        ).json())

except Exception as e:
    print(e, file=sys.stderr)
