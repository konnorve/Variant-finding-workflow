rule assembly_SPAdes:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
    output:
        scratch_dict["assembly"]["spades"] / "{sample}" / "scaffolds.fasta",
    conda:
        "../envs/spades.yaml"
    log:
        "logs/assemble_genomes/spades/{sample}.log"
    shell:
        "spades.py --threads {resources.tasks} -1 {input.r1} -2 {input.r2} -o $(dirname {output}) &> {log}"


rule assembly_metaFlye:
    input:
        lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read']
    output:
        scratch_dict["assembly"]["flye"] / "{sample}" / "assembly.fasta"
    conda:
        "../envs/flye.yaml"
    log:
        "logs/assemble_genomes/flye/{sample}.log"
    shell:
        "flye --threads {resources.tasks} --meta --pacbio-hifi {input} --out-dir $(dirname {output}) &> {log}"


def choose_assembly_method(wildcards):
    if wildcards.sample in PACBIO_SAMPLES:
        return scratch_dict["assembly"]["flye"] / f"{wildcards.sample}" / "assembly.fasta"
    elif wildcards.sample in ILLUMINA_SAMPLES:
        return scratch_dict["assembly"]["spades"] / f"{wildcards.sample}" / "scaffolds.fasta"
    else:
        raise ValueError("sample not in PacBio or Illumina Samples")
    

rule assembly_ragtag_scaffolding:
    input:
        assembly = choose_assembly_method,
        reference = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'genome_ref'],
    output:
        scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    conda:
        "../envs/ragtag.yaml"
    log:
        "logs/assemble_genomes/ragtag/{sample}.log"
    shell:
        "ragtag.py scaffold -w {input.reference} {input.assembly} -o $(dirname {output}) &> {log}"

