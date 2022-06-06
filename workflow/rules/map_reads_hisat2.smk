# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        index_dir = scratch_dict["genome_index"], 
        index = scratch_dict["done_files"]["index"]
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/hisat2/map_reads_PE/{sample}.log"
    shell:
        # HISAT2 Mapping:
        "hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output} 2> {log}"

rule make_sample_index_dir:
    output:
        directory(scratch_dict["genome_index"])
    shell:
        "mkdir -p {output}; sleep 1"

rule build_hisat2_index:
    input:
        out_dir = scratch_dict["genome_index"],
        ref = Path(config["input"]["genome_ref"]),
    output:
        touch(scratch_dict["done_files"]["index"])
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/hisat2/build_hisat2_index.log"
    shell:
        "hisat2-build {input.ref} {input.out_dir} > {log}"
        