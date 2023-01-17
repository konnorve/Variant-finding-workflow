
rule depth_histogram_data:
    input:
        depth = scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        genome_dup_stats = results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'depth_histogram_data', wildcards=[w.sample]),
        error = lambda w: mk_err('post_analysis', 'depth_histogram_data', wildcards=[w.sample]),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/gen_genome_depth_variant_data.py"

rule vcf_compare:
    input:
        expand(results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.{method}.vcf", sample=ALL_SAMPLES, method=config['variant_calling_methods']),
    output:
        directory(results_dict["vcf_compare"])
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'vcf_compare'),
        error = lambda w: mk_err('post_analysis', 'vcf_compare'),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/vcf_compare.py"

rule phased_gene_variants:
    input:
        vcf = results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        gff = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'gff_ref'],
    output:
        results_dict["phased_gene_variants"] / "{sample}.phased_gene_variants.tsv"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'phased_gene_variants', wildcards=[w.sample]),
        error = lambda w: mk_err('post_analysis', 'phased_gene_variants', wildcards=[w.sample]),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/get_phased_gene_variants.py"

rule lineage_occurance_table:
    input:
        phased_variants = expand(results_dict["phased_gene_variants"] / "{sample}.phased_gene_variants.tsv", sample=ALL_SAMPLES),
        sample_table = config['samples']
    output:
        results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_occurance.tsv"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'lineage_occurance_table'),
        error = lambda w: mk_err('post_analysis', 'lineage_occurance_table'),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/lineage_occurance.py"

rule heatmap:
    input:
        results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_occurance.tsv"
    output:
        png = results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_heatmap.png",
        html = results_dict["lineage_mutation_occurance"] / "lineage_gene_mutation_heatmap.html",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'heatmap'),
        error = lambda w: mk_err('post_analysis', 'heatmap'),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_heatmap.py"

rule depth_histogram_figure:
    input:
        results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    output:
        depth_histogram = results_dict["figures"]["depth_histogram"] / "{sample}_depth_histogram.png",
        depth_genome = results_dict["figures"]["depth_genome"] / "{sample}_depth_genome.png",
        depth_genome_zscore = results_dict["figures"]["depth_genome_zscore"] / "{sample}_depth_genome_zscore.png",
        coverage_gc = results_dict["figures"]["coverage_gc"] / "{sample}_coverage_gc.png",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'depth_histogram_figure', wildcards=[w.sample]),
        error = lambda w: mk_err('post_analysis', 'depth_histogram_figure', wildcards=[w.sample]),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_dup_del_plot.py"

rule dotplot:
    input:
        mums = results_dict["raw_data"]["mummer"] / "{sample}.mums" ,
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["figures"]["dotplot"] / "{sample}_dotplot.png"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('post_analysis', 'dotplot', wildcards=[w.sample]),
        error = lambda w: mk_err('post_analysis', 'dotplot', wildcards=[w.sample]),
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_dotplot.py"