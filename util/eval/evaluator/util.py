
_QUAST_TSV_TO_COLUMN_NAME = {
    '# contigs': 'n_contigs',
    'Largest contig': 'largest_contig',
    'N50': 'N50',
    'NG50': 'NG50',
    'NA50': 'NA50',
    'NGA50': 'NGA50',
    '# mismatches per 100 kbp': 'mismatches_per_100kbp',
    '# indels per 100 kbp': 'indels_per_100kbp',
    'Largest alignment': 'largest_alignment',
    'Total aligned length': 'total_aligned_length',
}


def parse_quast_tsv(file):
    ret = {}
    with open(file) as tsv:
        for line in tsv.readlines():
            tsv_name, data = line.strip().rsplit('\t')
            db_name = _QUAST_TSV_TO_COLUMN_NAME.get(tsv_name)
            if db_name is not None:
                ret[db_name] = data

    return ret
