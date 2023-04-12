rule call_variants_bcftools:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_standard/{sample}.log"
    shell:
        "bcftools mpileup --threads {resources.tasks} -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call --threads {resources.tasks} -mv -Ov -o {output.scratch} &> {log} && "
        "cp {output.scratch} {output.final}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_all/{sample}.log"
    shell:
        "bcftools mpileup --threads {resources.tasks} -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call --threads {resources.tasks} -cv -Ov -o {output.scratch} &> {log} && "
        "cp {output.scratch} {output.final}"


rule call_variants_freebayes:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/call_variants/freebayes/{sample}.log"
    shell:
        "freebayes -p 1 -f {input.ref_genome} {input.deduped_bam} > {output.scratch} 2> {log} && "
        "cp {output.scratch} {output.final}"
