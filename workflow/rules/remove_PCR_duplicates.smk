rule remove_PCR_duplicates:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        out_bam = scratch_dict["deduped"] / "{sample}_deduped.bam",
        dedup_metadata = scratch_dict["deduped"] / "{sample}_deduped_metadata.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates -I {input} -O {output.out_bam} -M {output.dedup_metadata} --REMOVE_SEQUENCING_DUPLICATES"
        