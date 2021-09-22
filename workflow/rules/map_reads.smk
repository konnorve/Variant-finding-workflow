# map sample reads to concatenated genome
rule map_reads:
    input:
        r1 = trimmed_reads_dir / "{sample}_1_trimmed.fastq",
        r2 = trimmed_reads_dir / "{sample}_2_trimmed.fastq",
        index_dir = genome_index_parent_dir / "{sample}" / "concat_genome", 
        index = done_file_dir / "build_hisat2_index.{sample}.done"
    output:
        temp(mapped_reads_dir / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    shell:
        # HISAT2 Mapping:
        """
        hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output}
        """

rule make_sample_index_dir:
    output:
        directory(genome_index_parent_dir / "{sample}" / "concat_genome")
    shell:
        "mkdir -p {output}"

rule build_hisat2_index:
    input:
        out_dir = genome_index_parent_dir / "{sample}" / "concat_genome",
        ref = reference_genome_file,
    output:
        touch(done_file_dir / "build_hisat2_index.{sample}.done")
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2-build {input.ref} {input.out_dir}"