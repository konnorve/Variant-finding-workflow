import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

depth_table = pd.read_table(snakemake.input[0], header=0, names=['contig', 'pos', 'depth'])

arr = depth_table.depth.to_numpy()

BIN_SIZE=10000

n = int(len(arr) / BIN_SIZE)

bin_size = int(len(arr) / n)

split_arr = []
count = 0
while count < len(arr):
    if count+bin_size > len(arr):
        split_arr.append(arr[count:])
    else:
        split_arr.append(arr[count:count+bin_size])
    count += bin_size

averaged = [np.mean(a) for a in split_arr]

plt.hist(averaged)

plt.savefig(snakemake.output["depth_histogram"])

plt.close()

plt.plot([(i*BIN_SIZE)+(BIN_SIZE/2) for i in range(len(averaged))], averaged)

plt.savefig(snakemake.output["depth_genome"])

plt.close()