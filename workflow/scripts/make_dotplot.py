from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Bio import SeqIO

# From mummer 3 manual: http://mummer.sourceforge.net/manual/
# the three columns list the position in the reference sequence, the position in the query sequence, and the length of the match respectively

def parse_mums_to_dict(mums_file, ref_contig=None):
    d = dict()
    with open(mums_file) as f:
        nl = next(f, None)
        while nl is not None:
            pl = nl
            nl = next(f, None)
            mums = []
            while nl is not None and ">" not in nl:
                split_line = nl.split()
                if len(split_line) == 3:
                    if ref_contig is None:
                        raise ValueError('ref_contig cannot be None')
                    else:
                        mum = {k:c(v) for k, c, v in zip(['ref_pos', 'query_pos', 'colinear_len'], [int]*3, split_line)}
                        mum['ref_contig'] = ref_contig
                        mums.append(mum)
                if len(split_line) == 4:
                    mum = {k:c(v) for k, c, v in zip(['ref_contig', 'ref_pos', 'query_pos', 'colinear_len'], [str, int, int, int], split_line)}
                    mums.append(mum)
                nl = next(f, None)
            if len(mums) > 0:
                d[pl.strip()[2:]] = mums
    return d

def plot_mums_plotly(mums_dict, ref_len_dict, query_len_dict, outpath):
    fig = go.Figure()

    query_ser = pd.Series(query_len_dict)
    query_ser = query_ser.cumsum() - query_ser

    ref_ser = pd.Series(ref_len_dict)
    ref_ser = ref_ser.iloc[::-1]
    ref_ser = ref_ser.cumsum() - ref_ser

    for i, (k, v) in enumerate(mums_dict.items()):
        x = []
        y = []
        for m in v:
            rc = m['ref_contig']
            rp = m['ref_pos'] + ref_ser.loc[rc]
            qp = m['query_pos'] + query_ser.loc[k]
            l = m['colinear_len']
            x.extend([rp, rp+l, None])
            y.extend([qp, qp+l, None])
        fig.add_trace(
            go.Scatter(x=x, y=y, name=k, mode='lines+markers')
        )

    fig.update_xaxes({
        'tickmode': 'array',
        'tickvals': list(ref_ser.values),
        'ticktext': list(ref_ser.index),
    })

    fig.update_yaxes({
        'tickmode': 'array',
        'tickvals': list(query_ser.values),
        'ticktext': list(query_ser.index),
    })
    
    fig.write_html(outpath)

def get_seqrecord_dict(fasta_path):
    with open(fasta_path) as in_seq_handle:
        return SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

def get_seq_length_dict(fasta_path):
    srd = get_seqrecord_dict(fasta_path)
    return {k:len(v.seq) for k,v in srd.items()}
    

mums_file = snakemake.input["mums"]
ref_genome = snakemake.input["ref_genome"]
query_genome = snakemake.input["query_genome"]
dotplot_outpath = snakemake.output[0]

ref_len_dict = get_seq_length_dict(ref_genome)
query_len_dict = get_seq_length_dict(query_genome)

if len(ref_len_dict) > 1:
    mums_dict = parse_mums_to_dict(mums_file)
elif len(ref_len_dict) == 1:
    mums_dict = parse_mums_to_dict(mums_file, ref_contig = list(ref_len_dict.keys())[0])
else:
    raise ValueError()

plot_mums_plotly(mums_dict, ref_len_dict, query_len_dict, dotplot_outpath)