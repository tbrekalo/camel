import sqlalchemy as sa
import sqlalchemy.orm as sa_orm


from .base import Base


class Quast(Base, ):
    __tablename__ = 'quast'

    id = sa.Column(sa.Integer, primary_key=True)
    n_contigs = sa.Column(sa.Integer, nullable=False)
    largest_contig = sa.Column(sa.Integer, nullable=False)
    N50 = sa.Column(sa.Integer, nullable=False)
    NG50 = sa.Column(sa.Integer, nullable=False)
    NA50 = sa.Column(sa.Integer, nullable=False)
    NGA50 = sa.Column(sa.Integer, nullable=False)
    mismatches_per_100kbp = sa.Column(sa.Float, nullable=False)
    indels_per_100kbp = sa.Column(sa.Float, nullable=False)
    largest_alignment = sa.Column(sa.Integer, nullable=False)
    total_aligned_length = sa.Column(sa.Integer, nullable=False)

    run_id = sa.Column(sa.Integer, sa.ForeignKey('run.id'))
    run = sa_orm.relationship('Run', back_populates='quast')
