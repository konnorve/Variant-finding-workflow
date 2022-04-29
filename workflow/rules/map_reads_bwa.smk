# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["index"]
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources: 
        mem_mb=100000,
    conda:
        "../envs/bwa.yaml"
    shell:
        # BWA Mapping:
        """
        bwa mem {input.ref} {input.r1} {input.r2} > {output.sam_out}
        """

rule index_genome:
    input:
        ref = Path(config["input"]["genome_ref"])
    output:
        touch(scratch_dict["done_files"]["index"]),
    resources:
        mem_mb=10000,
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.ref}"
