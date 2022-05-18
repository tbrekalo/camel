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


def plot_rng_region(pile: cp.Pile, read: cp.NucleicAcid,  dst_folder: Path, win_len: int = 250):
  fetch_fns: Dict[str, Callable[[cp.Coverage], np.uint16]] = {
      'A': lambda c: c.a,
      'C': lambda c: c.c,
      'G': lambda c: c.g,
      'T': lambda c: c.t,
      # 'mat': lambda c: c.match,
      # 'mis': lambda c: c.mismatch,
      'del': lambda c: c.deletion,
      'ins': lambda c: c.insertion,
  }

  slice_first = np.random.randint(0, len(pile.coverages) - win_len)
  slice_last = slice_first + win_len

  covgs: List[cp.Coverage] = pile.coverages[slice_first:slice_last]
  read_phred = np.frombuffer(read.inflate_quality().encode(), dtype=np.uint8)

  slice_phred = pd.DataFrame(
    data=read_phred[slice_first:slice_last],
    columns=['phred33'],
    index=np.arange(0, win_len))

  df = pd.DataFrame({'pos': np.arange(len(covgs))})
  for name, fetch_fn in fetch_fns.items():
    df[name] = np.fromiter(map(fetch_fn, covgs), dtype=np.uint16)
  df = df.melt(id_vars='pos', var_name='align_type', value_name='freq')

  fig, ax = plt.subplots(2, 1, figsize=(16, 12), dpi=500)

  plt_cov = sns.histplot(data=df, x='pos', hue='align_type', weights='freq',
                          binwidth=1, multiple='stack', ax=ax[0])
  plt_cov.set_ylabel('coverage') 

  plt_phred = sns.histplot(data=slice_phred, x=slice_phred.index, weights='phred33', 
    multiple='stack', binwidth=1, ax=ax[1])
  plt_phred.set_ylabel('phred33')
  plt_phred.axhline(np.median(read_phred), color='r', ls='--')

  dst_path = dst_folder.joinpath(f'pile_{pile.id:06d}')
  fig.savefig(dst_path)

  print(f'[camel::plot::plot_rng_region] saved {dst_path}', file=sys.stderr)

def plot_ins_distr(piles: List[cp.Pile], dst_folder: Path):
  ins_distr = np.array(cp.sample_ins_distr(piles))
  df = pd.DataFrame(data=ins_distr, columns=['freq', 'count'])[15:100]

  fig, ax = plt.subplots(figsize=(16,12), dpi=500)
  plt_ins_rate = sns.histplot(
    data=df, x='freq', weights='count',
    binwidth=1, ax=ax)

  plt_ins_rate.set_xlabel('insertion coverage value')
  plt_ins_rate.set_ylabel('count')

  dst_path = dst_folder.joinpath('ins_distr')
  fig.savefig(dst_path)

state = cp.init_state(32, './camel_log')
reads = cp.load_sequences(state, sys.argv[1:])
piles: List[cp.Pile] = cp.deserialize_piles(state, "./camel_piles")

print(f'[camel::plot] loaded {len(piles)} piles and {len(reads)} reads', file=sys.stderr)

read_index: Dict[str, int] = {}
for idx, read in enumerate(reads):
  read_index[read.name] = idx

dst_folder: Path = Path('./plots')
if dst_folder.exists():
  shutil.rmtree(dst_folder)
  dst_folder.mkdir(exist_ok=True)

for pile in random.choices(piles, k=10):
  plot_rng_region(pile, reads[read_index[pile.seq_name]], dst_folder,
                  win_len=min(250, int(len(pile.coverages) / 2)))
