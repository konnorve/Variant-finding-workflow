rule remove_PCR_duplicates:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        out_bam = scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
        dedup_metadata = scratch_dict["deduped"]["metadata"] / "{sample}_deduped_metadata.txt"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('remove_PCR_duplicates', 'remove_PCR_duplicates', wildcards=[w.sample]),
        error = lambda w: mk_err('remove_PCR_duplicates', 'remove_PCR_duplicates', wildcards=[w.sample]),
    conda:
        "../envs/picard.yaml"
    log: 
        "logs/remove_PCR_duplicates/picard/{sample}.log"
    shell:
        "picard MarkDuplicates -I {input} -O {output.out_bam} -M {output.dedup_metadata} --REMOVE_SEQUENCING_DUPLICATES &> {log}"
        

rule bam_coverage:
    input:
        scratch_dict["deduped"]["bams"] / "{sample}_deduped.bam",
    output:
        scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = lambda w: mk_out('remove_PCR_duplicates', 'bam_coverage', wildcards=[w.sample]),
        error = lambda w: mk_err('remove_PCR_duplicates', 'bam_coverage', wildcards=[w.sample]),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/remove_PCR_duplicates/samtools_depth/{sample}.log"
    shell:
        "samtools depth -o {output} {input}"