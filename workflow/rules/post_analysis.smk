

# rule vcf_annotater:
#     input:
#         vcf = scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf",
#         genbank = Path(config["input"]["genbank_ref"])
#     output:
#         results_dict["annotated_vcfs"] / "{sample}.{method}.vcf"
#     conda:
#         "../envs/vcf_annotater.yaml"
#     log: 
#         "logs/post_analysis/vcf_annotater/{sample}.{method}.log"
#     shell:
#         "vcf-annotator --output {output} {input.vcf} {input.genbank} > {log}"

rule depth_histogram_data:
    input:
        depth = scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv",
        ref = Path(config["input"]["genome_ref"]),
    output:
        genome_dup_stats = results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/gen_genome_depth_variant_data.py"

rule vcf_compare:
    input:
        expand(results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.{method}.vcf", sample=SAMPLES, method=config['variant_calling_methods']),
    output:
        directory(results_dict["vcf_compare"])
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/vcf_compare.py"

rule phased_gene_variants:
    input:
        vcf = results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        ref = Path(config["input"]["genome_ref"]),
        gff = Path(config["input"]["gff_ref"]),
    output:
        results_dict["phased_gene_variants"] / "{sample}.phased_gene_variants.tsv"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/get_phased_gene_variants.py"

rule depth_histogram_figure:
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

rule dotplot:
    input:
        mums = results_dict["raw_data"]["mummer"] / "{sample}.mums" 
        ref_genome = Path(config["input"]["genome_ref"]),
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["figures"]["dotplot"] / "{sample}_dotplot.png"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_dotplot.py"