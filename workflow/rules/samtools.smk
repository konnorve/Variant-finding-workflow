
rule convert_sam2bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/convert_sam2bam/{sample}.log"
    shell:
        "samtools view -S -b {input} > {output} 2> {log}"

rule sort_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam")
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort_bam/{sample}.log"
    shell:
        "samtools sort {input} -o {output} > {log}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index_bam/{sample}.log"
    shell:
        "samtools index -b {input} > {log}"

# rule index_fasta:
#     input:
#         reference_genome_file,
#     output:
#         reference_genome_index_file,
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools faidx {input}"
    