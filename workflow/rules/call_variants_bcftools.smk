rule call_variants_bcftools:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        variant_dir / "{sample}" / "{sample}_bcftools_standard.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -mv -Ov -o {output}"

rule call_variants_bcftools_all:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        variant_dir / "{sample}" / "{sample}_bcftools_all.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -c -Ov -o {output}"

rule link_bcftool_all:
    input:
        variant_dir / "{sample}" / "{sample}_bcftools_all.vcf"
    output:
        raw_all_calls_vcf_dir / "{sample}_bcftools_all.vcf"
    shell:
        "ln -s {input} {output}; sleep 1"

rule cut_metadata:
    input:
        raw_all_calls_vcf_dir / "{sample}_bcftools_all.vcf"
    output:
        all_calls_vcf_pre_df_dir / "{sample}_bcftools_all.tsv"
    shell:
        "grep -v '##' {input} > {output}; sleep 1"