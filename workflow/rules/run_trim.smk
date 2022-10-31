
rule run_trim:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'fwd_read'],
        r2 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'rev_read'],
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        o2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/run_trim/{sample}.log"
    shell:
        "bbduk.sh "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1 > {log}"