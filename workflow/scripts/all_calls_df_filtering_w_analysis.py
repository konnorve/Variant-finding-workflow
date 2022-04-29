from pathlib import Path
import deep_variant_calling as dvc

global_dir = Path("/nobackup1/kve/2021_Variant_Calling_Sean_Genomes")

reference_genome_file = global_dir / "input_data" / "reference_genome" / "SynechococcusRS9916.fna"
reference_genome_annotation_file = global_dir / "input_data" / "annotations" / "SynechococcusRS9916.gff"

variant_df_outpath = global_dir / "output_data" / 'all_calls_for_pandas' / 'variant_dfs' / 'all_genomes_all_variants.tsv'
filtered_variant_df_outpath = global_dir / "output_data" / 'all_calls_for_pandas' / 'variant_dfs' / 'all_genomes_filtered_variants.tsv'
variant_effect_df_outpath = global_dir / "output_data" / 'all_calls_for_pandas' / 'variant_dfs' / 'all_genomes_variant_effects.tsv'
alt_genome_dir = global_dir / "output_data" / "alt_genomes"

dvc.df_filtered_analysis(variant_df_outpath, reference_genome_annotation_file, reference_genome_file, filtered_variant_df_outpath, variant_effect_df_outpath, alt_genome_dir)