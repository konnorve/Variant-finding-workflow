from pathlib import Path
import deep_variant_calling as dvc
print("arr")
dvc.df_creation(snakemake.input, snakemake.output[0])