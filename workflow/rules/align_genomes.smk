rule align_genomes_mummer:
    input:
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
        query_genome = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    output:
        results_dict["raw_data"]["mummer"] / "{sample}.mums"
    conda:
        "../envs/mummer.yaml"
    log:
        "logs/align_genomes/mummer/{sample}.log"
    shell:
        "mummer {input.ref_genome} {input.query_genome} > {output} 2> {log}"