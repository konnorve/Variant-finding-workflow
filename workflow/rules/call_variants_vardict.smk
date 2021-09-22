rule call_variants_vardict:
    input:
        deduped_bam = deduped_dir / "{sample}_deduped.bam",
        deduped_index = mapped_reads_dir / "{sample}_mapped_sorted.bam.bai",
        ref_genome = reference_genome_file,
        ref_index = reference_genome_index_file
    output:
        vardict_dir / "{sample}_vardict_variants.vcf"
    resources:
        mem_mb=100000,
    conda:
        "../envs/vardict.yaml"
    params:
        allele_freq=0.05
    shell:
        "vardict -G {input.ref_genome} "
        "-f {params.allele_freq} -N {wildcards.sample} "
        "-b {input.deduped_bam} "
        "| teststrandbias.R | var2vcf_valid.pl "
        "-N {wildcards.sample} -E -f {params.allele_freq} > {output}"
        