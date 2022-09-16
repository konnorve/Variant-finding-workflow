

# rule vcf_annotater:
#     input:
#         vcf = scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf",
#         genbank = Path(config["input"]["genbank_ref"])
#     output:
#         results_dict["annotated_vcfs"] / "{sample}.{method}.vcf"
#     conda:
#         "../envs/vcf_annotater.yaml"
#     log: 
#         "logs/post_analysis/vcf_annotater/{sample}.{method}.log"
#     shell:
#         "vcf-annotator --output {output} {input.vcf} {input.genbank} > {log}"

rule depth_histogram:
    input:
        depth = scratch_dict["deduped"]["depth"] / "{sample}_depth.tsv",
        ref = Path(config["input"]["genome_ref"]),
    output:
        genome_dup_stats = results_dict["raw_data"]["genome_dup_stats"] / "{sample}_genome_bin_stats.tsv"
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/gen_genome_depth_variant_data.py"
