rule remove_PCR_duplicates:
    input:
        mapped_reads_dir / "{sample}_mapped_sorted.bam",
    output:
        out_bam = deduped_dir / "{sample}_deduped.bam",
        dedup_metadata = deduped_dir / "{sample}_deduped_metadata.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates -I {input} -O {output.out_bam} -M {output.dedup_metadata} --REMOVE_SEQUENCING_DUPLICATES"
        