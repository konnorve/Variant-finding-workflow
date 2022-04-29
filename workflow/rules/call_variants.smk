rule call_variants_bcftools:
    input:
        deduped_bam = scratch_dict["deduped"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_standard.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -mv -Ov -o {output}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = scratch_dict["deduped"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.bcftools_all.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -c -Ov -o {output}"


rule call_variants_freebayes:
    input:
        deduped_bam = scratch_dict["deduped"] / "{sample}_deduped.bam",
        ref_genome = Path(config["input"]["genome_ref"])
    output:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.freebayes.vcf"
    resources:
        mem_mb=32000,
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref_genome} {input.deduped_bam} > {output}"


rule link_bcftool_all:
    input:
        scratch_dict["vcfs"] / "{sample}" / "{sample}.{method}.vcf"
    output:
        scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf"
    shell:
        "ln -s {input} {output}; sleep 1"

rule cut_metadata:
    input:
        scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf"
    output:
        scratch_dict["all_calls"]["calls_only"] / "{sample}.{method}.tsv"
    shell:
        "grep -v '##' {input} > {output}; sleep 1"
