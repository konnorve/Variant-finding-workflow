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
