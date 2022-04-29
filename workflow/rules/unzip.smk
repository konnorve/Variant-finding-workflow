rule unzip:
    input:
        scratch_dict["trimmed_reads"] / "{anything}.fastq.gz",
    output:
        scratch_dict["trimmed_reads"] / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"
        