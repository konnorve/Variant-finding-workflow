
rule map_reads_PacBio:
    input:
        r = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read'],
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["index"]
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources: 
        mem_mb=100000,
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/bwa/map_reads_PE/{sample}.log"
    shell:
        # BWA Mapping:
        "bwa mem {input.ref} {input.r} > {output.sam_out} 2> {log}"

rule index_genome:
    input:
        Path(config["input"]["genome_ref"])
    output:
        scratch_dict["done_files"]["index"]
    resources:
        mem_mb=10000,
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/bwa/index_genome.log"
    shell:
        "bwa index {input} &> {log}; touch {output}; sleep 1"
