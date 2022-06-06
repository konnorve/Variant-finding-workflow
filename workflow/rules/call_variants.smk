rule call_variants_bcftools:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_standard/{sample}.log"
    shell:
        "bcftools mpileup -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -mv -Ov -o {output} &> {log}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/call_variants/bcftools_all/{sample}.log"
    shell:
        "bcftools mpileup -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -c -Ov -o {output} &> {log}"


rule call_variants_freebayes:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.freebayes.vcf"
    resources:
        mem_mb=32000,
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/call_variants/freebayes/{sample}.log"
    shell:
        "freebayes -f {input.ref_genome} {input.deduped_bam} > {output} 2> {log}"


rule link_bcftool_all:
    input:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.{method}.vcf"
    output:
        scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf"
    log:
        "logs/call_variants/link_bcftool_all/{sample}.{method}.log"
    shell:
        "ln -s {input} {output}; sleep 1"

rule cut_metadata:
    input:
        scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf"
    output:
        scratch_dict["all_calls"]["calls_only"] / "{sample}.{method}.tsv"
    shell:
        "grep -v '##' {input} > {output}; sleep 1"
