rule de_novo_SPAdes:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
    output:
        outdir = directory(scratch_dict["assembly"]["spades"] / "{sample}"),
        assembly = scratch_dict["assembly"]["spades"] / "{sample}" / "scaffolds.fasta",
    resources:
        mem_mb=100000,
    conda:
        "../envs/spades.yaml"
    log:
        "logs/assemble_genomes/spades/{sample}.log"
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o {output.outdir} &> {log}"

rule ragtag_scaffolding:
    input:
        assembly = scratch_dict["assembly"]["spades"] / "{sample}" / "scaffolds.fasta",
        reference = Path(config["input"]["genome_ref"]),
    output:
        outdir = directory(scratch_dict["assembly"]["ragtag"] / "{sample}"),
        assembly = scratch_dict["assembly"]["ragtag"] / "{sample}" / "ragtag.scaffold.fasta",
    resources:
        mem_mb=100000,
    conda:
        "../envs/ragtag.yaml"
    log:
        "logs/assemble_genomes/ragtag/{sample}.log"
    shell:
        "ragtag.py scaffold -w {input.reference} {input.assembly} -o {output.outdir} &> {log}"

rule nucmer_align:
    input:
        assembly = scratch_dict["assembly"]["spades"] / "{sample}" / "scaffolds.fasta",
        reference = Path(config["input"]["genome_ref"]),
    output:
        scratch_dict["assembly"]["ragtag"] / "{sample}.delta"
    conda:
        "../envs/mummer.yaml"
    params:
        prefix=str(scratch_dict["assembly"]["ragtag"]) + "/{sample}"
    log: 
        "logs/assemble_genomes/nucmer_align/{sample}.log"
    shell:
        "nucmer --prefix={params.prefix} {input.assembly} {input.reference} &> {log}"

rule nucmer_coords:
    input:
        scratch_dict["assembly"]["ragtag"] / "{sample}.delta"
    output:
        scratch_dict["assembly"]["ragtag"] / "{sample}.coords"
    conda:
        "../envs/mummer.yaml"
    log: 
        "logs/assemble_genomes/nucmer_coords/{sample}.log"
    shell:
        "show-coords -r -c -l {input} > {output} 2> {log}"

