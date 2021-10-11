from pathlib import Path
import deep_variant_calling as dvc

io_dir = Path('/nobackup1/kve/2021_Variant_Calling_Sean_Genomes/output_data/all_calls_for_pandas')
all_calls_dir = io_dir / 'vcf_pre_df'
variant_dfs_dir = io_dir / 'variant_dfs'
variant_df_outpath = variant_dfs_dir / 'all_genomes_all_variants.tsv'

dvc.df_creation(all_calls_dir, variant_df_outpath)
