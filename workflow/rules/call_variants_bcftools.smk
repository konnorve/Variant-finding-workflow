rule call_variants_bcftools:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        bcftools_dir / "{sample}_bcftools_variants.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -mv -Ov -o {output}"

rule call_variants_bcftools_felix:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        bcftools_felix_dir / "{sample}_bcftools_felix_variants.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools mpileup -Ou -q30 -x -d3000 -f {input.ref_genome} {input.deduped_bam} | "
        "bcftools call -c -Ov -o {output}"