
rule post_analaysis_depth_histogram_data:
    input:
        depth = scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        genome_dup_stats = results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/gen_genome_depth_variant_data.py"

rule post_analaysis_vcf_compare:
    input:
        expand(results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.{method}.vcf", sample=ALL_SAMPLES, method=config['variant_calling_methods']),
    output:
        directory(results_dict["vcf_compare"])
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/vcf_compare.py"

rule post_analaysis_phased_gene_variants:
    input:
        vcf = results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        gff = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'gff_ref'],
    output:
        results_dict["phased_gene_variants"] / "{sample}.phased_gene_variants.tsv"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/get_phased_gene_variants.py"

rule post_analaysis_lineage_occurance_table:
    input:
        phased_variants = expand(results_dict["phased_gene_variants"] / "{sample}.phased_gene_variants.tsv", sample=ALL_SAMPLES),
        sample_table = config['samples']
    output:
        results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_occurance.tsv"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/lineage_occurance.py"

rule post_analaysis_heatmap:
    input:
        results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_occurance.tsv"
    output:
        png = results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_heatmap.png",
        html = results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_heatmap.html",
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_heatmap.py"

rule post_analaysis_depth_histogram_figure:
    input:
        results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    output:
        depth_histogram = results_dict["figures"]["depth_histogram"] / "{sample}_depth_histogram.png",
        depth_genome = results_dict["figures"]["depth_genome"] / "{sample}_depth_genome.png",
        depth_genome_zscore = results_dict["figures"]["depth_genome_zscore"] / "{sample}_depth_genome_zscore.png",
        coverage_gc = results_dict["figures"]["coverage_gc"] / "{sample}_coverage_gc.png",
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_dup_del_plot.py"

rule post_analaysis_dotplot:
    input:
        mums = results_dict["raw_data"]["mummer"] / "{sample}.mums" ,
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["figures"]["dotplot"] / "{sample}_dotplot.html"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_dotplot.py"
