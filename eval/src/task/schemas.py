import re

from pydantic import BaseModel, FilePath, validator

from task.exceptions import (TaskConfigInvalidErrorThreshold,
                             TaskConfigInvalidExeError, TaskRunInvalidValue)
from utils import to_camel


class BaseConfig(BaseModel):
    class Config:
        allow_mutation = False
        alias_generator = to_camel
        allow_population_by_field_name = True


class TaskConfig(BaseConfig):
    exe: str
    threads: int
    window_length: int
    error_threshold: float
    reads_path: FilePath
    overlaps_path: FilePath

    @validator('exe')
    @classmethod
    def exe_validator(cls, value):
        if not re.match(r'\w+', value):
            raise TaskConfigInvalidExeError(
                value,
                "invalid exe value"
            )
        return value

    @validator('error_threshold')
    @classmethod
    def error_threshold_validator(cls, value):
        if value < 0 or value >= 1:
            raise TaskConfigInvalidErrorThreshold(
                value,
                "error_threshold should be between 0 and 1"
            )

        return value


class TaskRun(BaseConfig):
    peak_memory_mib: float
    runtime_s: float

    @validator('peak_memory_mib')
    @classmethod
    def peak_memory_mib_validator(cls, value):
        if value < 0:
            raise TaskRunInvalidValue(
                value,
                "peak_memory_mib has to be greater or equal to 0"
            )

        return round(value, 2)

    @validator('runtime_s')
    @classmethod
    def runtime_s_validator(cls, value):
        if value < 0:
            raise TaskRunInvalidValue(
                value,
                "runtime_s has to be greater or equal to 0"
            )

        return round(value, 2)

    class Config:
        allow_mutation = False
        alias_generator = to_camel
        allow_population_by_field_name = True


class TaskInfo(BaseConfig):
    task_config: TaskConfig
    task_run: TaskRun

    class Config:
        allow_mutation = False
        alias_generator = to_camel
        allow_population_by_field_name = True
