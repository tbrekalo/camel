from argparse import ArgumentParser
from typing import Any, Dict, List, Optional


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
        self._add_unary_arg('--database', type=str, default=['pyeval.db'])
        self._add_unary_arg('--work-dir', type=str, default=['/tmp/pyeval'])


class _RaconArgBuilder(_ArgBuilder):
    def __init__(self):
        super().__init__('correction')
        self._add_unary_arg('--threads', type=int, default=1)
        self._add_unary_arg('--window-length', type=int, default=500)
        self._add_unary_arg('--error-threshold', type=float, default=0.3)
        self._argparser.add_argument('overlaps')
        self._argparser.add_argument('reads')


class _CamelArgBuilder(_RaconArgBuilder):
    def __init__(self):
        super().__init__()

        self._add_unary_arg('--n-overlaps', type=int, default=128)
        self._add_unary_arg('--comment', type=str, required=True)


def parse_args(args: Optional[List[str]] = None) -> Dict[str, Any]:
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

    exe_args = exe_args_parser.parse_args(eval_args.args.split(' '))
    return {
        'executable': eval_args.executable,
        'work_dir': eval_args.work_dir,
        'db_path': eval_args.database,
    } | vars(exe_args)
