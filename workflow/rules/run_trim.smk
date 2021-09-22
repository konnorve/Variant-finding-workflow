
rule run_trim:
    input:
        r1 = resequenced_genome_reads_dir / "{sample}1.fastq.gz",
        r2 = resequenced_genome_reads_dir / "{sample}2.fastq.gz",
        ref = adapter_file,
    output:
        o1 = trimmed_reads_dir / "{sample}_1_trimmed.fastq.gz",
        o2 = trimmed_reads_dir / "{sample}_2_trimmed.fastq.gz",
    resources:
        mem_mb=100000,
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"