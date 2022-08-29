from pathlib import Path
import deep_variant_calling as dvc
import shutil
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

gff_path = snakemake.input["gff_path"]
ref_path = snakemake.input["ref_path"]
full_variant_df_path = snakemake.input["vcf_df_path"]

output_dir = Path(snakemake.output["outdir"])
method = output_dir.name

filtered_variant_df_outpath = output_dir / f"all_genomes_filtered_variants.{method}.tsv"
intragenic_figure_dir = output_dir / "intragenic" / "figures"
intragenic_df_outpath = output_dir / "intragenic" / f"all_genomes_intragenic_variants.{method}.tsv"
alt_genome_dir = output_dir / "intragenic" / "alt_genomes"
intergenic_figure_dir = output_dir / "intergenic" / "figures"
intergenic_df_outpath = output_dir / "intergenic" / f"all_genomes_intergenic_variants.{method}.tsv"

for d in [intragenic_figure_dir, alt_genome_dir, intergenic_figure_dir]:
    d.mkdir(parents=True, exist_ok=True)

full_variant_df = dvc.readVariantDF(full_variant_df_path)

if method == "bcftools_all":
    filtered_variant_df = dvc.filterOutNonSigVariants(full_variant_df)
    dvc.saveVariantDF(filtered_variant_df, filtered_variant_df_outpath)
else:
    dvc.saveVariantDF(full_variant_df, filtered_variant_df_outpath)
    filtered_variant_df = full_variant_df

# filtered
dvc.intragenic_analysis(filtered_variant_df, gff_path, ref_path, intragenic_figure_dir, intragenic_df_outpath, alt_genome_dir)
dvc.intergenic_analysis(filtered_variant_df, gff_path, ref_path, intergenic_figure_dir, intergenic_df_outpath)

gene_obj_type = type(snakemake.output["genomes"])
logging.info(f"type genomes obj: {gene_obj_type}")
for f in snakemake.output["genomes"]:
    if not Path(f).exists():
        logging.info(f)
        shutil.copyfile(ref_path, f)