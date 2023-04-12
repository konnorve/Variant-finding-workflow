# map sample reads to concatenated genome
rule mapping_bowtie2_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "illumina_bowtie2" / "{sample}" / "{sample}_mapped.sam"),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {input.ref} $(dirname {output})/index
        bowtie2 -p {resources.tasks} -x $(dirname {output})/index -1 {input.r1} -2 {input.r2} -S {output.sam_out}
        rm $(dirname {output})/index*
        """

rule mapping_bwa_PacBio:
    input:
        r = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read'],
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "pacbio_bwa" / "{sample}" / "{sample}_mapped.sam"),
    conda:
        "../envs/bwa.yaml"
    shell:
        # BWA Mapping:
        """
        bwa index -p $(dirname {output})/index {input.ref}
        bwa mem -t {resources.tasks} $(dirname {output})/index {input.r} > {output.sam_out}
        rm $(dirname {output})/index*
        """

# map sample reads to concatenated genome
rule mapping_bwa_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "illumina_bwa" / "{sample}" / "{sample}_mapped.sam"),
    conda:
        "../envs/bwa.yaml"
    shell:
        # BWA Mapping:
        """
        bwa index -p $(dirname {output})/index {input.ref}
        bwa mem -t {resources.tasks} $(dirname {output})/index {input.r1} {input.r2} > {output.sam_out}
        rm $(dirname {output})/index*
        """
