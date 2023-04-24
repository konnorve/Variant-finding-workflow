

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
        var_pos = scratch_dict["variant positions"] / "{reference}.vcf"
    run:
        # combines all VCFs into one file
        vcfs = []
        for vcf in input:
            temp_df = pd.read_table(vcf, comment="#", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "other"])
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
        df = df.sort_values(['CHROM','POS','sample']).set_index(['CHROM','POS','sample'])
        # add filtering steps here for positions
        # saves positions to table
        index = df.index.droplevel("sample").drop_duplicates()
        index.to_frame().to_csv(output[0], index=False, header=False, sep='\t')


def get_var_pos_file(wildcards):
    reference = SAMPLE_TABLE.loc[wildcards.sample, 'reference']
    return scratch_dict["variant positions"] / f"{reference}.vcf"

rule call_variants_samtools_pileup:
    input:
        bam = scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
        ref_genome = get_reference,
        var_pos = get_var_pos_file,
    output:
        mpileup=scratch_dict["mpileup"] / "{sample}.pileup"
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
        "--positions {input.var_pos} "
        "--fasta-ref {input.ref_genome} "
        "{input.bam} >> {output.mpileup}"

def get_reference_mpileups(wildcards):
    reference_samples = SAMPLE_TABLE[SAMPLE_TABLE['reference']==wildcards.reference].index.to_list()
    return expand(scratch_dict["mpileup"] / "{sample}.pileup", sample=reference_samples)

rule filter_samtools_mpileup:
    input:
        pileups = get_reference_mpileups,
    output:
        table = scratch_dict["filtered mpileups"] / "{reference}.tsv"
    run:
        # combines all pileups into one file
        pileups = []
        for pileup in input.pileups:
            temp_df = pd.read_table(pileup)
            temp_df['sample'] = Path(pileup).stem
            pileups.append(temp_df)
        df = pd.concat(pileups)
        # indexes dataframe on scaffold, position, and sample
        df = df.sort_values(['chrom','pos','sample']).set_index(['chrom','pos','sample'])
        # saves positions to table
        df.to_csv(output.table, sep='\t')
