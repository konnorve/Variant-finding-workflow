import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqUtils
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

df = pd.read_table(snakemake.input["???"])

plt.hist(df.mean_bin_depth)
plt.savefig(snakemake.output["depth_histogram"])
plt.close()

plt.scatter(df[df['2sigma']==True].bin_midpoint, df[df['2sigma']==True].normalized_corrected_bin_depth_change, c='r', marker='.')
plt.scatter(df[df['2sigma']==False].bin_midpoint, df[df['2sigma']==False].normalized_corrected_bin_depth_change, c='b', marker='.')
plt.savefig(snakemake.output["depth_genome"])
plt.close()

plt.scatter(df[df['2sigma']==True].bin_midpoint, df[df['2sigma']==True].normalized_corrected_bin_depth, c='r', marker='.')
plt.scatter(df[df['2sigma']==False].bin_midpoint, df[df['2sigma']==False].normalized_corrected_bin_depth, c='b', marker='.')
plt.savefig(snakemake.output["depth_genome_zscore"])
plt.close()


plt.scatter(df.bin_gc, df.normalized_bin_depth, c='b', marker='.')
plt.plot(df.bin_gc, [m*x+b for x in df.bin_gc], color='red')
plt.savefig(snakemake.output["coverage_gc"])
plt.close()