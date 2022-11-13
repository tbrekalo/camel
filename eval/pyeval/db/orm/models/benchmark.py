import sqlalchemy as sa
import sqlalchemy.orm as sa_orm

from .base import Base


class Benchmark(Base,):
    __tablename__ = 'benchmark'

    id = sa.Column(sa.Integer, primary_key=True)
    threads = sa.Column(sa.Integer, nullable=False)
    runtime_s = sa.Column(sa.Integer, nullable=False)
    peak_memory_mib = sa.Column(sa.Integer, nullable=False)

    run_id = sa.Column(sa.Integer, sa.ForeignKey('run.id'))
    run = sa_orm.relationship('Run', back_populates='benchmark')
