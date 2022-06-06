# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["index"]
    output:
        sam_out = temp(output_path_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources: 
        mem_mb=100000,
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/bowtie2/map_reads_PE/{sample}.log"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output.sam_out} &> {log}"

rule index_genome:
    input:
        ref = Path(config["input"]["genome_ref"])
    output:
        touch(scratch_dict["done_files"]["index"]),
    resources:
        mem_mb=10000,
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/bowtie2/index_genome.log"
    shell:
        "bowtie2-build {input.ref} {input.ref} &> {log}"

