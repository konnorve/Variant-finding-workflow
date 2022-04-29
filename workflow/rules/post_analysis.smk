
rule generate_all_variant_calls_df:
    input:
        expand(scratch_dict["all_calls"]["calls_only"] / "{sample}.{method}.tsv", sample=SAMPLES, method=config['variant_calling_methods'])
    output:
        results_dict["all_variants_df"]
    conda:
        "../envs/blah.yaml"
    script:
        "../scripts/generate_counts_metadata_dfs.py"


rule vcf_annotater:
    input:
        vcf = scratch_dict["all_calls"]["raw"] / "{sample}.{method}.vcf",
        genbank = Path(config["input"]["genbank_ref"])
    output:
        results_dict["annotated_vcfs"] / "{sample}.{method}.vcf"
    conda:
        "../envs/vcf_annotater.yaml"
    shell:
        "vcf-annotator.py --output {output} {input.vcf} {input.genbank}"
