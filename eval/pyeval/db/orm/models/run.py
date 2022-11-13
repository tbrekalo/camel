import sqlalchemy as sa
import sqlalchemy.orm as sa_orm

from .base import Base


class Run(Base, ):
    __tablename__ = 'run'

    id = sa.Column(sa.Integer, primary_key=True)
    executable = sa.Column(sa.String)

    args = sa_orm.relationship('Args', back_populates='run', uselist=False)
    benchmark = sa_orm.relationship(
        'Benchmark', back_populates='run', uselist=False)
    quast = sa_orm.relationship('Quast', back_populates='run', uselist=False)
