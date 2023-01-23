import asyncio
import sys

from args import parse_args
from tasks import eval_correction
from db.context import create_engine, DBContext


async def main():
    eval_cfg = parse_args()
    engine = create_engine(eval_cfg.db_path)
    with DBContext(engine) as context:
        await eval_correction(context, eval_cfg)

try:
    asyncio.run(main())
except Exception as err:
    print(f'[pyeval] {err}', file=sys.stderr)
