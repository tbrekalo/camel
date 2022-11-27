import sqlalchemy as sa
import sqlalchemy.orm as sa_orm

from .base import Base


class Args(Base, ):
    __tablename__ = 'args'

    id = sa.Column(sa.Integer, primary_key=True)
    window_length = sa.Column(sa.Integer, nullable=False)
    error_threshold = sa.Column(sa.Float, nullable=False)
    reads = sa.Column(sa.String, nullable=False)
    overlaps = sa.Column(sa.String, nullable=False)
    n_overlaps = sa.Column(sa.Integer, nullable=True)

    run_id = sa.Column(sa.Integer, sa.ForeignKey('run.id'))
    run = sa_orm.relationship('Run', back_populates='args')
