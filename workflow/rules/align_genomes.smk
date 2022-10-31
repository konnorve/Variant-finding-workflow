rule run_mummer:
    input:
        ref_genome = Path(config["input"]["genome_ref"]),
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["raw_data"]["mummer"] / "{sample}.mums" 
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '1G',
        ntasks = 1,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/mummer.yaml"
    log:
        "logs/align_genomes/mummer/{sample}.log"
    shell:
        "mummer {input.ref_genome} {input.query_genome} > {output} 2> {log}"