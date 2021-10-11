rule analyze_variants_vep:
    input:
        vcf = variant_dir / "{sample}" / "{sample}_{variant_call_tool}.vcf",
        annotation = new_annotation_file,
        annotation_index = new_annotation_index,
        reference = reference_genome_file

    output:
        variant_dir / "{sample}" / "{sample}_{variant_call_tool}_vep_results.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/ensembl_vep.yaml"
    shell:
        "vep -i {input.vcf} -o {output} --sift b --gff {input.annotation} --fasta {input.reference}" 

rule prep_gff_vep:
    input:
        reference_genome_annotation_file
    output:
        new_annotation_file
    resources:
        mem_mb=100000,
    conda:
        "../envs/ensembl_vep.yaml"
    shell:
        'grep -v "#" {input} | sort -k1,1 -k4,4n -k5,5n -t$"\t" | bgzip -c > {output}'

rule index_gff:
    input:
        new_annotation_file
    output:
        new_annotation_index
    resources:
        mem_mb=100000,
    conda:
        "../envs/ensembl_vep.yaml"
    shell:
        'tabix -p gff {input}'