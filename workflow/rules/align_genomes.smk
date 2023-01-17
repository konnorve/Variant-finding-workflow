rule run_mummer:
    input:
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["raw_data"]["mummer"] / "{sample}.mums" 
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '1G',
        ntasks = 1,
        time = '1-0',
        output = lambda w: mk_out('align_genomes', 'run_mummer', wildcards=[w.sample]),
        error = lambda w: mk_err('align_genomes', 'run_mummer', wildcards=[w.sample]),
    conda:
        "../envs/mummer.yaml"
    log:
        "logs/align_genomes/mummer/{sample}.log"
    shell:
        "mummer {input.ref_genome} {input.query_genome} > {output} 2> {log}"