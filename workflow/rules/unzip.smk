rule unzip:
    input:
        trimmed_reads_dir / "{anything}.fastq.gz",
    output:
        trimmed_reads_dir / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"
        