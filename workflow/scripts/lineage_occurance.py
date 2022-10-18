from pathlib import Path
import pandas as pd

def main(phased_variant_paths, sample_table_path, gene_proportions_outpath):
    sample_df = pd.read_table(sample_table_path, index_col='sample_name')

    phased_variant_paths = pd.concat({Path(f).name.split('.')[0]: pd.read_table(f, index_col=['seq_id', 'ID']) for f in phased_variant_paths}, names=['sample','chrom','gene'])
    phased_variant_paths = phased_variant_paths[phased_variant_paths['synonymous']==False]

    variants = phased_variant_paths.index.to_frame(index=False)
    variants = variants.join(sample_df[['phenotype','lineage']], on='sample')
    # variants['lineage'] = variants['sample'].map(sample_df['lineage'].to_dict())

    lineage_variants = variants.groupby(['chrom','gene','phenotype','lineage']).size()
    lineage_variants = lineage_variants >= 2
    lineage_variants = lineage_variants[lineage_variants]

    treatment_variants = lineage_variants.index.to_frame(index=False)
    treatment_variants = treatment_variants.groupby(['chrom','gene','phenotype']).size()

    unique_phenotype_lineages = dict()
    for k, v in set(zip(sample_df['phenotype'], sample_df['lineage'])):
        unique_phenotype_lineages[k] = unique_phenotype_lineages.setdefault(k, 0) + 1

    phenotype_variants = treatment_variants.unstack('phenotype').fillna(0)
    phenotype_variants = phenotype_variants.div(list(map(unique_phenotype_lineages.get, phenotype_variants.columns)), axis='columns')

    phenotype_variants.reset_index().to_csv(gene_proportions_outpath, sep='\t', index=False)

phased_variant_paths = list(Path('/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/results/phased_gene_variants').glob("*.tsv"))
sample_table_path = "/home/kve/scripts/variant-calling-workflow/config/samples.tsv"

main(phased_variant_paths, sample_table_path)
# main(snakemake.input["vcf"], snakemake.input["ref"], snakemake.input["gff"], snakemake.output[0])