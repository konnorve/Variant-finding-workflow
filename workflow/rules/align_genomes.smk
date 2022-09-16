rule run_mummer:
    input:
        ref_genome = Path(config["input"]["genome_ref"]),
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["raw_data"]["mummer"] / "{sample}.mums" 
    resources:
        mem_mb=100000,
    conda:
        "../envs/mummer.yaml"
    log:
        "logs/align_genomes/mummer/{sample}.log"
    shell:
        "mummer {input.ref_genome} {input.query_genome} > {output} 2> {log}"