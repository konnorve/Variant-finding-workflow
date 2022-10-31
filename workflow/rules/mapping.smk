# map sample reads to concatenated genome
rule map_reads_bowtie2_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["bowtie2_index"]
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "illumina_bowtie2" / "{sample}_mapped.sam"),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/mapping/map_reads_bowtie2_PE/{sample}.log"
    shell:
        "bowtie2 -p {resources.ntasks} -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output.sam_out} &> {log}"

rule index_genome_bowtie2:
    input:
        ref = Path(config["input"]["genome_ref"])
    output:
        touch(scratch_dict["done_files"]["bowtie2_index"]),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/mapping/index_genome.log"
    shell:
        "bowtie2-build {input.ref} {input.ref} &> {log}"

rule map_reads_bwa_PacBio:
    input:
        r = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'ccs_read'],
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["bwa_index"]
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "pacbio_bwa" / "{sample}_mapped.sam"),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/mapping/map_reads_bwa_pacbio/{sample}.log"
    shell:
        # BWA Mapping:
        "bwa mem -t {resources.ntasks} {input.ref} {input.r} > {output.sam_out} 2> {log}"

# map sample reads to concatenated genome
rule map_reads_bwa_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = Path(config["input"]["genome_ref"]),
        indexing = scratch_dict["done_files"]["bwa_index"]
    output:
        sam_out = temp(scratch_dict["mapped_reads"] / "illumina_bwa" / "{sample}_mapped.sam"),
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/mapping/map_reads_bwa_PE/{sample}.log"
    shell:
        # BWA Mapping:
        "bwa mem -t {resources.ntasks} {input.ref} {input.r1} {input.r2} > {output.sam_out} 2> {log}"

rule index_genome_bwa:
    input:
        Path(config["input"]["genome_ref"])
    output:
        touch(scratch_dict["done_files"]["bwa_index"])
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/bwa.yaml"
    log: 
        "logs/mapping/index_genome.log"
    shell:
        "bwa index {input} &> {log}; sleep 10"

# map sample reads to concatenated genome
rule map_reads_hisat2_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        index_dir = scratch_dict["genome_index"], 
        index = scratch_dict["done_files"]["hisat_index"]
    output:
        temp(scratch_dict["mapped_reads"] / "illumina_hisat2" / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/mapping/map_reads_hisat2_PE/{sample}.log"
    shell:
        # HISAT2 Mapping:
        "hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output} 2> {log}"

rule make_hisat2_index_dir:
    output:
        directory(scratch_dict["genome_index"])
    shell:
        "mkdir -p {output}; sleep 1"

rule build_hisat2_index:
    input:
        out_dir = scratch_dict["genome_index"],
        ref = Path(config["input"]["genome_ref"]),
    output:
        touch(scratch_dict["done_files"]["hisat_index"])
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/mapping/index_genome.log"
    shell:
        "hisat2-build {input.ref} {input.out_dir} &> {log}"
        