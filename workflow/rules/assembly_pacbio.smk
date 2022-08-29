rule metaFlye:
    input:
        lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read']
    output:
        outdir = directory(scratch_dict["assembly"]["flye"] / "{sample}"),
        assembly = scratch_dict["assembly"]["flye"] / "{sample}" / "assembly.fasta",
    resources:
        mem_mb=100000,
    conda:
        "../envs/flye.yaml"
    log:
        "logs/assemble_genomes/flye/{sample}.log"
    shell:
        "flye --meta --pacbio-hifi {input} --out-dir {output.outdir} &> {log}"


rule ragtag_scaffolding:
    input:
        assembly = scratch_dict["assembly"]["flye"] / "{sample}" / "assembly.fasta",
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
        assembly = scratch_dict["assembly"]["flye"] / "{sample}" / "assembly.fasta",
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

