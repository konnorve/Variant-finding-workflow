
def choose_mapped_inputs(wildcards):
    if wildcards.sample in PACBIO_SAMPLES:
        method = "pacbio_bwa"
    elif wildcards.sample in ILLUMINA_SAMPLES:
        method = f"illumina_{config['mapper']}"
    else:
        raise ValueError("sample not in PacBio or Illumina Samples")
    return scratch_dict["mapped_reads"] / method / "{sample}_mapped.sam",


rule convert_sam2bam:
    input:
        choose_mapped_inputs
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
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
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort_bam/{sample}.log"
    shell:
        "samtools sort {input} -o {output} &> {log}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index_bam/{sample}.log"
    shell:
        "samtools index -b {input} &> {log}"

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
    