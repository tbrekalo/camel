import argparse

from pathlib import Path

import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='plot window len distribution' )
  parser.add_argument('--src_dir', type=str, nargs=1)
  parser.add_argument('--dst_dir', type=str, nargs=1)

  args = parser.parse_args()
  src_dir = Path(args.src_dir[0])
  dst_dir = Path(args.dst_dir[0])
  
  def parse_content(p: Path):
    with open(p) as src:
      line = src.readline().strip()
      if len(line) != 0:
        return np.fromiter(
          map(int, line.split()), dtype=np.uint32)
      return []

  win_len_logs = [Path(e) for e in src_dir.glob('window_lens_*')]
  win_lens, cnts =  np.unique(np.fromiter(
    (x for p in win_len_logs for x in parse_content(p)), 
      dtype=np.uint32), return_counts=True)

  df = pd.DataFrame(data=np.array([win_lens, cnts]).T,
    columns=['win_len', 'freq'])

  print(df)

  fig, ax = plt.subplots(figsize=(16,12), dpi=500)
  sns.histplot(data=df, x='win_len', weights='freq', binwidth=1, ax=ax)

  dst_path = dst_dir.joinpath('win_len_distr')
  fig.savefig(dst_path)
