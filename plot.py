import camelpy as cp

from typing import List

tp_handle = cp.create_thread_pool(1)
piles: List[cp.Pile] = cp.deserialize_piles(tp_handle, "./test_dump")
