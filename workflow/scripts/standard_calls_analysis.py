from pathlib import Path
import deep_variant_calling as dvc

global_dir = Path("/nobackup1/kve/2021_Variant_Calling_Sean_Genomes")

input_dir = global_dir / "input_data"
output_dir = global_dir / "output_data"
results_dir = global_dir / "results"

resequenced_genome_reads_dir = input_dir / "resequenced_genome_reads"
reference_genome_file = input_dir / "reference_genome" / "SynechococcusRS9916.fna"
reference_genome_index_file  = input_dir / "reference_genome" / "SynechococcusRS9916.fna.fai"
adapter_file = input_dir / "adapters" / "all_illumina_adapters.fa"
reference_genome_annotation_file = input_dir / "annotations" / "SynechococcusRS9916.gff"

io_dir = output_dir / 'all_calls_for_pandas'
all_calls_dir = io_dir / 'vcf_pre_df'
variant_dfs_dir = io_dir / 'variant_dfs'
filtered_variant_df_path = variant_dfs_dir / 'standard_calls_filtered.tsv'
variant_effect_df_outpath = variant_dfs_dir / 'all_genomes_standard_variant_effects.tsv'
alt_genome_dir = output_dir / "alt_genomes"

dvc.df_standard_analysis(filtered_variant_df_path, reference_genome_annotation_file, reference_genome_file, variant_effect_df_outpath, alt_genome_dir)