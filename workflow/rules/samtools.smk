
def choose_mapped_inputs(wildcards):
    if wildcards.sample in PACBIO_SAMPLES:
        method = "pacbio_bwa"
    elif wildcards.sample in ILLUMINA_SAMPLES:
        method = f"illumina_{config['mapper']}"
    else:
        raise ValueError("sample not in PacBio or Illumina Samples")
    return scratch_dict["mapped_reads"] / method / "{sample}" / "{sample}_mapped.sam",


rule convert_sam2bam:
    input:
        choose_mapped_inputs
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input}"
