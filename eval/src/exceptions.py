class ReferenceNotFound(Exception):
    def __init__(self, value: str, message: str) -> None:
        super().__init__(message)
        self.value = value
