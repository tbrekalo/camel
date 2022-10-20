from argparse import ArgumentParser
from typing import Any, Dict, List, Optional
from pathlib import Path

from config import EvalConfig


class _ArgBuilder:
    def __init__(self, name: Optional[str] = None):
        if name is None:
            self._argparser = ArgumentParser()
        else:
            self._argparser = ArgumentParser(name)

    def _add_unary_arg(self, *args, **kwargs) -> None:
        self._argparser.add_argument(
            *args, nargs='?', action='store', **kwargs)

    def build(self):
        return self._argparser


class _EvalArgBuilder(_ArgBuilder):
    def __init__(self):
        super().__init__(name='pyeval')

        self._add_unary_arg('--executable', type=str, required=True)
        self._add_unary_arg('--args', type=str, required=True)
        self._add_unary_arg('--threads', type=int, default=1)
        self._add_unary_arg('--database', type=str, default='pyeval.db')
        self._add_unary_arg('--work-dir', type=str, default='/tmp/pyeval')
        self._add_unary_arg('--reference', type=str, required=True)
        self._add_unary_arg('--comment', type=str)


class _RaconArgBuilder(_ArgBuilder):
    def __init__(self):
        super().__init__('correction')
        self._add_unary_arg('--window-length', type=int, default=500)
        self._add_unary_arg('--error-threshold', type=float, default=0.3)
        self._argparser.add_argument('reads')
        self._argparser.add_argument('overlaps')


class _CamelArgBuilder(_RaconArgBuilder):
    def __init__(self):
        super().__init__()

        self._add_unary_arg('--n-overlaps', type=int, default=128)


def parse_args(args: Optional[List[str]] = None) -> EvalConfig:
    eval_args_parser = _EvalArgBuilder().build()

    if args:
        print(args)
        eval_args = eval_args_parser.parse_args(args)
    else:
        eval_args = eval_args_parser.parse_args()

    if eval_args.executable == 'camel':
        exe_args_parser = _CamelArgBuilder().build()
    elif eval_args.executable == 'racon':
        exe_args_parser = _RaconArgBuilder().build()
    else:
        raise RuntimeError('[pyeval] unsupported')

    exe_args = {'threads': eval_args.threads}
    parsed_args = vars(exe_args_parser.parse_args(eval_args.args.split(' ')))
    for k, v in parsed_args.items():
        exe_args[f'{k.replace("_", "-")}'] = v

    return EvalConfig(
        Path(eval_args.work_dir),
        eval_args.executable,
        exe_args,
        eval_args.threads,
        Path(eval_args.reference),
        Path(eval_args.database),
        eval_args.comment)
