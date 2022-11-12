import asyncio
import sys

import db

from args import parse_args
from job import eval_correction


async def main():
    eval_cfg = parse_args()
    quast_data, runtime_s, peak_memory_MiB = await eval_correction(eval_cfg)

    con = db.get_db_connection(eval_cfg.db_path)
    db.create_table_for_exe(
        con, eval_cfg.executable, eval_cfg.args)

    data = {
        'runtime_s': int(runtime_s),
        'peak_memory_MiB': peak_memory_MiB,
    } | {
        k: str(v) for k, v in eval_cfg.args.items()
    } | quast_data

    db.insert_dict_to_runs_table(con, eval_cfg.executable, data)


try:
    asyncio.run(main())
except Exception as err:
    print(f'[pyeval] {err}', file=sys.stderr)
