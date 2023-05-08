rule call_variants_samtools_pileup:
    input:
        bam = scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
        ref_genome = get_reference,
    output:
        mpileup=temp(scratch_dict["mpileup"] / "{sample}.pileup")
    conda:
        "../envs/samtools.yaml"
    params:
        max_depth=config['secondary variant calling']['max depth'],
        mapping_qual_thresh=config['secondary variant calling']['mapping quality threshold'],
        base_qual_thresh=config['secondary variant calling']['base quality threshold'],
    shell:
        "echo -e 'chrom\tpos\tref\tcov\tbases\tquals\tMAPQ\ttail_dist' > {output.mpileup} && "
        "samtools mpileup "
        "--output-MQ --output-BP "
        "--min-MQ {params.mapping_qual_thresh} "
        "--min-BQ {params.base_qual_thresh} "
        "--max-depth {params.max_depth} "
        "--fasta-ref {input.ref_genome} "
        "--no-output-ins --no-output-ins "
        "--no-output-del --no-output-del "
        "--no-output-ends "
        "{input.bam} >> {output.mpileup}"

rule gzip_mpileup:
    input:
        scratch_dict["mpileup"] / "{sample}.pileup"
    output:
        scratch_dict["mpileup"] / "{sample}.pileup.gz"
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip {input}"
    

rule call_variants_bcftools_pileup:
    input:
        bam = scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
        ref_genome = get_reference,
    output:
        vcf = scratch_dict["vcfs"] / "{sample}.vcf"
    conda:
        "../envs/bcftools.yaml"
    params:
        max_depth=config['initial variant calling']['max depth'],
        mapping_qual_thresh=config['initial variant calling']['mapping quality threshold'],
        base_qual_thresh=config['initial variant calling']['base quality threshold'],
        snp2indeldist=config['initial variant calling']['min distance of snp to indel'],
    shell:
        "bcftools mpileup -Ou "
        "-q {params.mapping_qual_thresh} "
        "-Q {params.base_qual_thresh} "
        "-d {params.max_depth} "
        "-f {input.ref_genome} "
        "{input.bam} | "
        "bcftools call -mv -Ou | "
        "bcftools filter --SnpGap {params.snp2indeldist} > {output.vcf}"
        
        # "--full-BAQ "

def get_reference_vcfs(wildcards):
    reference_samples = SAMPLE_TABLE[SAMPLE_TABLE['reference']==wildcards.reference].index.to_list()
    return expand(scratch_dict["vcfs"] / "{sample}.vcf", sample=reference_samples)

rule variant_union:
    input:
        get_reference_vcfs
    output:
        var_pos = scratch_dict["variant positions"] / "{reference}.pos"
    run:
        # combines all VCFs into one file
        vcfs = []
        for vcf in input:
            temp_df = pd.read_table(vcf, comment="#", names=["chrom", "pos", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "other"])
            temp_df['sample'] = Path(vcf).stem
            vcfs.append(temp_df)
        df = pd.concat(vcfs)
        print(df)
        # filters out ambiguous positions
        df = df[~df['ALT'].str.contains(",")]
        # only keeps SNPs
        df = df[df['REF'].apply(lambda s: len(s)==1)]
        df = df[df['ALT'].apply(lambda s: len(s)==1)]
        # indexes dataframe on scaffold, position, and sample
        df = df.sort_values(['chrom','pos','sample']).set_index(['chrom','pos','sample'])
        # add filtering steps here for positions
        # saves positions to table
        index = df.index.droplevel("sample").drop_duplicates()
        index.to_frame().to_csv(output[0], index=False, sep='\t')

def get_reference_mpileups(wildcards):
    reference_samples = SAMPLE_TABLE[SAMPLE_TABLE['reference']==wildcards.reference].index.to_list()
    return expand(scratch_dict["mpileup"] / "{sample}.pileup.gz", sample=reference_samples)

rule make_candidate_mutation_table:
    input:
        pileups = get_reference_mpileups,
        var_pos = scratch_dict["variant positions"] / "{reference}.pos"
    output:
        table = scratch_dict["filtered mpileups"] / "{reference}.tsv.gz"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/make_candidate_mutation_table.py"

rule decorate_candidate_mutation_table:
    input:
        scratch_dict["filtered mpileups"] / "{reference}.tsv.gz"
    output:
        results_dir / "{reference}_candidate_mutation_table.tsv.gz"
    conda:
        "../envs/post_analysis.yaml"
    script:
        "../scripts/decorate_candidate_mutation_table.py"