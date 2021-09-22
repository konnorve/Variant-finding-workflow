
rule convert_sam2bam:
    input:
        mapped_reads_dir / "{sample}_mapped.sam",
    output:
        temp(mapped_reads_dir / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        mapped_reads_dir / "{sample}_mapped.bam",
    output:
        mapped_reads_dir / "{sample}_mapped_sorted.bam",  # TODO: once we can confirm workflow we should change this to temp
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule index_bam:
    input:
        mapped_reads_dir / "{sample}_mapped_sorted.bam",
    output:
        temp(mapped_reads_dir / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input}"

rule index_fasta:
    input:
        reference_genome_file,
    output:
        reference_genome_index_file,
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"
    