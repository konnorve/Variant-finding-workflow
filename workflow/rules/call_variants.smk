rule call_variants_bcftools:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf",
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_standard/{sample}.log"
    shell:
        "bcftools mpileup -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -mv -Ov -o {output.scratch} &> {log} && "
        "ln {output.scratch} {output.final}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf",
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_all/{sample}.log"
    shell:
        "bcftools mpileup -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -cv -Ov -o {output.scratch} &> {log} && "
        "ln {output.scratch} {output.final}"


rule call_variants_freebayes:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.freebayes.vcf",
    resources:
        mem_mb=32000,
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/call_variants/freebayes/{sample}.log"
    shell:
        "freebayes -p 1 -f {input.ref_genome} {input.deduped_bam} > {output.scratch} 2> {log} && "
        "ln {output.scratch} {output.final}"

rule call_variants_breseq_consensus:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        outdir=directory(scratch_dict["vcfs"] / "{sample}" / "{sample}.breseq_consensus"),
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.breseq_consensus" / "output" / "output.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.breseq_consensus.vcf",
    resources:
        mem_mb=120000,
    threads: 10
    conda:
        "../envs/breseq.yaml"
    log:
        "logs/call_variants/breseq_consensus/{sample}.log"
    shell:
        "breseq -j {threads} -r {input.ref_genome} -o {output.outdir} {input.r1} {input.r1} &> {log} && "
        "ln {output.scratch} {output.final}"

rule call_variants_breseq_polymorphism:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        outdir=directory(scratch_dict["vcfs"] / "{sample}" / "{sample}.breseq_polymorphism"),
        scratch=scratch_dict["vcfs"] / "{sample}" / "{sample}.breseq_polymorphism" / "output" / "output.vcf",
        final=results_dict["raw_data"]["raw_vcfs"] / "{sample}" / "{sample}.breseq_polymorphism.vcf",
    resources:
        mem_mb=120000,
    threads: 10
    conda:
        "../envs/breseq.yaml"
    log:
        "logs/call_variants/breseq_polymorphism/{sample}.log"
    shell:
        "breseq -p -j {threads} -r {input.ref_genome} -o {output.outdir} {input.r1} {input.r1} &> {log} && "
        "ln {output.scratch} {output.final}"
