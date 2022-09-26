import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqUtils
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

depth_table = pd.read_table(snakemake.input["depth"], header=0, names=['contig', 'pos', 'depth'])

with open(snakemake.input["ref"]) as in_seq_handle:
    sequence = next(SeqIO.parse(in_seq_handle, "fasta")).seq

depth_arr = depth_table.depth.to_numpy()

BIN_SIZE=1000

GENOME_SIZE=len(depth_arr)

bins = np.linspace(0, GENOME_SIZE, int(GENOME_SIZE/BIN_SIZE), dtype=int)

mean_genome_depth = np.mean(depth_arr)

bin_info = []
for i in range(len(bins)-1):
    bin_start = bins[i]
    bin_end = bins[i+1]
    bin_size = bin_end - bin_start
    bin_midpoint = (bin_end + bin_start) / 2
    mean_bin_depth = np.mean(depth_arr[bin_start:bin_end])
    normalized_bin_depth = mean_bin_depth / mean_genome_depth
    bin_seq = sequence[bin_start:bin_end]
    bin_gc = SeqUtils.GC(bin_seq)

    bin_info.append([bin_start, bin_end, bin_size, bin_midpoint, mean_bin_depth, normalized_bin_depth, bin_gc])

df = pd.DataFrame(bin_info, columns=["bin_start", "bin_end", "bin_size", "bin_midpoint", "mean_bin_depth", "normalized_bin_depth", "bin_gc"])

m, b = np.polyfit(df.bin_gc, df.normalized_bin_depth, 1)

df['gc_regression_prediction'] = df['bin_gc'].apply(lambda x: m*x+b)
df['normalized_corrected_bin_depth'] = df['normalized_bin_depth'] - df['gc_regression_prediction']
df['normalized_corrected_bin_depth_change'] = df['normalized_bin_depth'] / df['gc_regression_prediction']

genome_depth_sd = np.std(df.normalized_corrected_bin_depth)
df['z_score'] = df['normalized_corrected_bin_depth'] / genome_depth_sd
df['2sigma'] = df['z_score'].apply(lambda x: x > 2 or x < -2)

df.to_csv(snakemake.output["genome_dup_stats"], sep='\t')