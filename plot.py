import random
import shutil
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from typing import List

from camelpy import create_thread_pool, deserialize_piles, Pile


def plot_as_df(pile: Pile, win_len: int = 250) -> None:
  value_fetch_fns = {
      'match': lambda cov: cov.match,
      'deletion': lambda cov: cov.deletion,
      'insertion': lambda cov: cov.insertion,
      'mismatch': lambda cov: cov.mismatch
  }

  slice_first: int = np.random.randint(len(pile.coverages) - win_len)
  slice_last: int = slice_first + win_len

  covg_slice = pile.coverages[slice_first:slice_last]

  df = pd.DataFrame({'pos': np.arange(0, len(covg_slice))})
  for col_name, fetch_fn in value_fetch_fns.items():
    df.insert(0, col_name,
              np.fromiter(map(fetch_fn, covg_slice), dtype=np.uint16))

  df = df.melt('pos', var_name='alignment type', value_name='freq')
  print(df)

  plt.figure(figsize=(12.8, 9.6), dpi=500)
  plt.suptitle(f'({pile.id:06d} {pile.seq_name})')
  sns.histplot(data=df, x='pos', y='freq', stat='count', 
    hue='alignment type', palette='pastel',
    binwidth=1, multiple='stack', element='step')
  plt.savefig(f'./plots/{pile.id:06d}.png', format='png')


tp_handle = create_thread_pool(32)
piles: List[Pile] = deserialize_piles(tp_handle, './data')
print(f'loaded {len(piles)}')

if os.path.exists('./plots'):
  shutil.rmtree('./plots')
os.mkdir('./plots')

sns.set_theme('notebook')
sns.set_style('dark')

for pile in random.choices(piles, k=5):
  if len(pile.coverages) < 2000:
    continue
  plot_as_df(pile)
