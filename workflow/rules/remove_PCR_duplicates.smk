rule remove_PCR_duplicates:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        out_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        dedup_metadata = scratch_dict["deduped"]["metadata"] / "{sample}_deduped_metadata.txt"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates -I {input} -O {output.out_bam} -M {output.dedup_metadata} --REMOVE_SEQUENCING_DUPLICATES"
        

rule bam_coverage:
    input:
        scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
    output:
        scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools depth -o {output} {input}"