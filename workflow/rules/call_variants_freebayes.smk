rule call_variants_freebayes:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        variant_dir / "{sample}" / "{sample}_freebayes.vcf"
    resources:
        mem_mb=32000,
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref_genome} {input.deduped_bam} > {output}"
        