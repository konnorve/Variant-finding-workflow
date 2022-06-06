import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import seaborn as sns; sns.set()

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

logging.info(f"finished loading packages")

variant_df = pd.read_table(snakemake.input[0])

variant_df["pos_ref_alt"] = variant_df[["POS", "REF", "ALT"]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

unique_variants = sorted(list(variant_df["pos_ref_alt"].unique()))
unique_genomes = sorted(list(variant_df["genome_name"].unique()))

logging.info(f"unique_variants:\n{unique_variants}")
logging.info(f"unique_genomes:\n{unique_genomes}")

arr = np.zeros((len(unique_variants), len(unique_genomes)))

for g, v in zip(variant_df["genome_name"], variant_df["pos_ref_alt"]):
    arr[unique_variants.index(v), unique_genomes.index(g)] = 1

distance_matrix = pd.DataFrame(arr, columns=unique_genomes, index=unique_variants)

distance_matrix.to_csv(snakemake.output["matrix"], sep="\t")

logging.info(f"distance_matrix:\n{distance_matrix}")

if len(unique_variants)/3 > 50:
    figure = sns.clustermap(distance_matrix, figsize=(max([len(t) for t in unique_variants])/5+len(unique_genomes)/3, 50), col_cluster=True, row_cluster=True)
else:
    figure = sns.clustermap(distance_matrix, figsize=(max([len(t) for t in unique_variants])/5+len(unique_genomes)/2, len(unique_variants)/3), col_cluster=True, row_cluster=True)


ax = plt.gcf().axes[2]
text_list = ax.get_yticklabels()
ax.set_yticklabels([l.get_text()[:25] for l in text_list])

logging.info(f"sum of array: {np.sum(arr)}")
logging.info(f"len of variant df: {len(variant_df)}")

figure.savefig(snakemake.output["clustering"])
