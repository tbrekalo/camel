from pathlib import Path

import sqlalchemy as sa
import sqlalchemy.orm as sa_orm

from .orm.models.base import Base

from .orm.models.args import Args
from .orm.models.benchmark import Benchmark
from .orm.models.quast import Quast
from .orm.models.run import Run

def create_uri(sqlite_path: Path) -> str:
    return f'sqlite:///{sqlite_path}'

def create_engine(sqlite_path: Path) -> sa.engine.Engine:
    return sa.create_engine(f'sqlite:///{sqlite_path}', echo=False)


class DBContext:
    def __init__(self, engine: sa.engine.Engine):
        self.__engine = engine

        self.executable: str
        self.version: str
        self.comment: str

        self.args: Args
        self.benchmark: Benchmark
        self.quast: Quast

    def __enter__(self):
        Session = sa_orm.sessionmaker(bind=self.__engine)
        self.__session = Session()

        Base.metadata.create_all(self.__engine)
        return self

    def __exit__(self, err_type, value, traceback):
        if err_type is None:
            run = Run(
                executable=self.executable,
                version=self.version,
                comment=self.comment)
            self.__session.add(run)
            self.__session.commit()

            self.args.run_id = run.id
            self.benchmark.run_id = run.id
            self.quast.run_id = run.id

            self.__session.add(self.args)
            self.__session.add(self.benchmark)
            self.__session.add(self.quast)
            self.__session.commit()
        else:
            return True
