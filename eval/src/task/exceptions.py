class TaskConfigInvalidErrorThreshold(Exception):
    def __init__(self, value: str, message: str) -> None:
        super().__init__(message)
        self.value = value


class TaskConfigInvalidExeError(Exception):
    def __init__(self, value: str, message: str) -> None:
        super().__init__(message)
        self.value = value


class TaskRunInvalidValue(Exception):
    def __init__(self, value: str, message: str) -> None:
        super().__init__(message)
        self.value = value
