rule de_novo_SPAdes:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
    output:
        scratch_dict["assembly"]["spades"] / "{sample}" / "scaffolds.fasta",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('assembly', 'de_novo_SPAdes', wildcards=[w.sample]),
        error = lambda w: mk_err('assembly', 'de_novo_SPAdes', wildcards=[w.sample]),
    conda:
        "../envs/spades.yaml"
    log:
        "logs/assemble_genomes/spades/{sample}.log"
    shell:
        "spades.py --threads {resources.ntasks} -1 {input.r1} -2 {input.r2} -o $(dirname {output}) &> {log}"


rule metaFlye:
    input:
        lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read']
    output:
        scratch_dict["assembly"]["flye"] / "{sample}" / "assembly.fasta",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('assembly', 'metaFlye', wildcards=[w.sample]),
        error = lambda w: mk_err('assembly', 'metaFlye', wildcards=[w.sample]),
    conda:
        "../envs/flye.yaml"
    log:
        "logs/assemble_genomes/flye/{sample}.log"
    shell:
        "flye --threads {resources.ntasks} --meta --pacbio-hifi {input} --out-dir $(dirname {output}) &> {log}"


def choose_assembly_method(wildcards):
    if wildcards.sample in PACBIO_SAMPLES:
        return scratch_dict["assembly"]["flye"] / f"{wildcards.sample}" / "assembly.fasta"
    elif wildcards.sample in ILLUMINA_SAMPLES:
        return scratch_dict["assembly"]["spades"] / f"{wildcards.sample}" / "scaffolds.fasta"
    else:
        raise ValueError("sample not in PacBio or Illumina Samples")
    

rule ragtag_scaffolding:
    input:
        assembly = choose_assembly_method,
        reference = Path(config["input"]["genome_ref"]),
    output:
        scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('assembly', 'ragtag_scaffolding', wildcards=[w.sample]),
        error = lambda w: mk_err('assembly', 'ragtag_scaffolding', wildcards=[w.sample]),
    conda:
        "../envs/ragtag.yaml"
    log:
        "logs/assemble_genomes/ragtag/{sample}.log"
    shell:
        "ragtag.py scaffold -w {input.reference} {input.assembly} -o $(dirname {output}) &> {log}"

