rule call_variants_freebayes:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        ref_genome = reference_genome_file
    output:
        freebayes_dir / "{sample}_freebayes_variants.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref_genome} {input.deduped_bam} > {output}"
        