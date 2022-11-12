from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Union

StrongDict = Dict[str, Union[int, float, str]]


@dataclass
class EvalConfig:
    work_dir: Path

    executable: str
    args: StrongDict
    threads: int

    reference_path: Path
    db_path: Path

    comment: str
