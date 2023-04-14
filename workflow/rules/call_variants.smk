
rule call_variants_samtools_pileup:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        mpileup=scratch_dict["mpileup"] / "{sample}.samtools.pileup"
    conda:
        "../envs/samtools.yaml"
    params:
        max_depth=config['mpileup']['max depth'],
        mapping_qual_thresh=config['mpileup']['mapping quality threshold'],
        base_qual_thresh=config['mpileup']['base quality threshold'],
    shell:
        "echo -e 'chrom\tpos\tref\tcov\tbases\tquals\tMAPQ\ttail_dist'> {output.mpileup} && "
        "samtools mpileup "
        "--output-MQ --output-BP "
        "-q {params.mapping_qual_thresh} "
        "-Q {params.base_qual_thresh} "
        "-d {params.max_depth} -f {input.ref_genome} {input.deduped_bam} >> {output.mpileup}"



rule call_variants_bcftools_pileup:
    input:
        deduped_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        ref_genome = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        vcf=scratch_dict["mpileup"] / "{sample}.bcftools.pileup"
    conda:
        "../envs/bcftools.yaml"
    params:
        max_depth=config['mpileup']['max depth'],
        mapping_qual_thresh=config['mpileup']['mapping quality threshold'],
        base_qual_thresh=config['mpileup']['base quality threshold'],
    shell:
        "bcftools mpileup "
        "-q {params.mapping_qual_thresh} "
        "-Q {params.base_qual_thresh} "
        "-d {params.max_depth} "
        "-f {input.ref_genome} "
        "--output {output.vcf} --output-type v "
        "{input.deduped_bam}"
        
        # "--full-BAQ "
