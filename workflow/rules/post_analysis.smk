
rule generate_all_variant_calls_df:
    input:
        expand(scratch_dict["all_calls"]["calls_only"] / "{sample}.{{method}}.tsv", sample=SAMPLES)
    output:
        results_dict["all_variants_df"] / "all_genomes_all_variants.{method}.tsv"
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/all_calls_df_creation.py"

rule vcf_annotater:
    input:
        vcf = scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf",
        genbank = Path(config["input"]["genbank_ref"])
    output:
        results_dict["annotated_vcfs"] / "{sample}.{method}.vcf"
    conda:
        "../envs/vcf_annotater.yaml"
    log: 
        "logs/post_analysis/vcf_annotater/{sample}.{method}.log"
    shell:
        "vcf-annotator --output {output} {input.vcf} {input.genbank} > {log}"


rule variant_effects:
    input:
        ref_path = Path(config["input"]["genome_ref"]),
        gff_path = Path(config["input"]["gff_ref"]),
        vcf_df_path = results_dict["all_variants_df"] / "all_genomes_all_variants.{method}.tsv",
    output:
        outdir = directory(results_dir / "{method}"),
        genomes = expand(results_dir / "{{method}}" / "intragenic" / "alt_genomes" / "{sample}_alt.fasta", sample=SAMPLES),
        filtered_variant_df_outpath = results_dir / "{method}" / "all_genomes_filtered_variants.{method}.tsv",
        done = touch(scratch_dict["done_files"]["analysis"] / "{method}.done")
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/variant_effects.py"

rule depth_histogram:
    input:
        depth = scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv",
        ref = Path(config["input"]["genome_ref"]),
    output:
        depth_histogram = results_dict["depth_histograms"] / "{sample}_depth_histogram.png",
        depth_genome = results_dict["depth_genome"] / "{sample}_depth_genome.png",
        depth_genome_zscore = results_dict["depth_genome_zscore"] / "{sample}_depth_genome_zscore.png",
        genome_bin_stats = results_dict["genome_bin_stats"] / "{sample}_genome_bin_stats.tsv",
        coverage_gc = results_dict["coverage_gc_correlation"] / "{sample}_coverage_gc_corr.png"
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/depth_histogram.py"

rule make_dotplot:
    input:
        mums = scratch_dict["pairwise_alignment"] / "{sample}" / "mummer.mums",
        ref_genome = Path(config["input"]["genome_ref"]),
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["dotplots"] / "{sample}_dotplot.html"
    conda:
        "../envs/deep_variant_calling.yaml"
    log: 
        "logs/post_analysis/make_dotplot/{sample}.log"
    script:
        "../scripts/make_dotplot.py"


rule distance_metrics:
    input:
        results_dir / "{method}" / "all_genomes_filtered_variants.{method}.tsv"
    output:
        matrix = results_dict["distance_matix"] / "{method}_distance_matrix.tsv",
        clustering = results_dict["distance_clustering"] / "{method}_distance_clustering.png",
    conda:
        "../envs/deep_variant_calling.yaml"
    log: 
        "logs/post_analysis/distance_metrics/{method}.log"
    script:
        "../scripts/distance_metrics.py"
