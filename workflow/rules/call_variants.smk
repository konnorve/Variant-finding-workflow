rule call_variants_bcftools:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('call_variants', 'call_variants_bcftools', wildcards=[w.sample]),
        error = lambda w: mk_err('call_variants', 'call_variants_bcftools', wildcards=[w.sample]),
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_standard/{sample}.log"
    shell:
        "bcftools mpileup --threads {resources.ntasks} -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call --threads {resources.ntasks} -mv -Ov -o {output.scratch} &> {log} && "
        "cp {output.scratch} {output.final}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('call_variants', 'call_variants_bcftools_all', wildcards=[w.sample]),
        error = lambda w: mk_err('call_variants', 'call_variants_bcftools_all', wildcards=[w.sample]),
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_all/{sample}.log"
    shell:
        "bcftools mpileup --threads {resources.ntasks} -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call --threads {resources.ntasks} -cv -Ov -o {output.scratch} &> {log} && "
        "cp {output.scratch} {output.final}"


rule call_variants_freebayes:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('call_variants', 'call_variants_freebayes', wildcards=[w.sample]),
        error = lambda w: mk_err('call_variants', 'call_variants_freebayes', wildcards=[w.sample]),
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/call_variants/freebayes/{sample}.log"
    shell:
        "freebayes -p 1 -f {input.ref_genome} {input.deduped_bam} > {output.scratch} 2> {log} && "
        "cp {output.scratch} {output.final}"
