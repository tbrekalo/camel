import random
import shutil
import sys

import numpy as np
import pandas as pd
import camelpy as cp
import seaborn as sns

from pathlib import Path
from matplotlib import pyplot as plt
from typing import List, Dict, Callable


def plot_rng_region(pile: cp.Pile, dst_folder: Path, win_len: int = 250):
  fetch_fns: Dict[str, Callable[[cp.Coverage], np.uint16]] = {
      'match': lambda c: c.match,
      'deletion': lambda c: c.deletion,
      'insertion': lambda c: c.insertion,
      'mismatch': lambda c: c.mismatch
  }

  slice_first = np.random.randint(0, len(pile.coverages) - win_len)
  slice_last = slice_first + win_len

  covgs: List[cp.Coverage] = pile.coverages[slice_first:slice_last]

  df = pd.DataFrame({'pos': np.arange(len(covgs))})
  for name, fetch_fn in fetch_fns.items():
    df[name] = np.fromiter(map(fetch_fn, covgs), dtype=np.uint16)
  df = df.melt(id_vars='pos', var_name='align_type', value_name='freq')

  fig, ax = plt.subplots(figsize=(16, 12), dpi=500)
  sns.histplot(data=df, x='pos', hue='align_type', weights='freq',
                          binwidth=1, multiple='stack', ax=ax)

  dst_path = dst_folder.joinpath(f'pile_{pile.id:06d}')
  fig.savefig(dst_path)

  print(f'[camel::plot::plot_rng_region] saved {dst_path}', file=sys.stderr)


tp_handle = cp.create_thread_pool(32)
piles: List[cp.Pile] = cp.deserialize_piles(tp_handle, "./data")

print(f'[camel::plot] loaded {len(piles)} piles', file=sys.stderr)

dst_folder: Path = Path('./plots')
if dst_folder.exists():
  shutil.rmtree(dst_folder)
  dst_folder.mkdir(exist_ok=True)

for pile in random.choices(piles, k=3):
  plot_rng_region(pile, dst_folder,
                  win_len=min(250, int(len(pile.coverages) / 2)))
