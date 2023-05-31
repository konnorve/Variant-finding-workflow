# Snakemake workflow: Variant Analysis

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Konnor von Emster (@kve)

## Purpose

Used for detecting variants of multiple types.

1. Long (5-10kb) insertions, deletions, and transversions
2. Short (500nt-5kb) duplications, deletions
3. Variant effect on genes

## Outputs

* Long (5-10kb) insertions, deletions, and transversions
  - Dot/mummer plot
* Short (500nt-5kb) duplications, deletions
  - Scatter plot and linear regression of bin GC content vs Coverage
  - Scatter plot of bin genomic position vs normalized coverage
  - Scatter plot of bin genomic position vs coverage Z score
* Variant effect on genes
  - Table of phased variant effect on genes
  - Plots of each variant within it's genomic context
* Meta outputs
  - Phased variant occurance between control and treatment
  - ??? #TODO: Improve

## Process

General processes are as follows:

For each sample:

1. Trim adapters and poor quality base calls from reads
2. Map reads to genome
   1. Get coverage information (used for short duplications and deletions)
   2. Call short indels and snps using variant calling tool of choice (bcftools, freebayes, etc.)
      1. Applies variants to genome
      2. Applies variants to GFF
      3. Compares differences from gene to gene (misense, nonsense, frameshift, silent)
      4. Reports phased gene variants
3. Assembly of reads into genomes
   1. Mummer of assembled reads and genome
   2. Creates Interactive Dot Plot using Plotly

Meta processes involving all samples:
Phased gene variant occurance

## Rule Details

* align_genomes.smk
    - align_genomes_mummer
      - Uses Mummer to create alignments
* assembly.smk
    - assembly_SPAdes
      - Used for assembly of Illumina Reads
    - assembly_metaFlye
      - Used for assembly of PacBio Reads
    - assembly_ragtag_scaffolding
      - Used for scaffolding of final assemblies #TODO: add specificity to this description
* call_variants.smk
    - call_variants_bcftools
      - Variant calling method (creates VCF file) that uses bcftool's standard base calling methods
    - call_variants_bcftools_all
      - Variant calling method (creates VCF file) that uses bcftool's inclusive base calling method
    - call_variants_freebayes
      - Variant calling method (creates VCF file) that uses Freebayes
* mapping.smk
    - mapping_bowtie2_PE
        - Used for mapping Illumina Paired End reads onto reference genome using Bowtie2
    - mapping_bwa_PacBio
        - Used for mapping PacBio reads onto reference genome using BWA
    - mapping_bwa_PE
        - Used for mapping Illumina Paired End reads onto reference genome using BWA
* post_analysis.smk
    - post_analaysis_depth_histogram_data
        - bins mapping coverage data and creates statistics
    - post_analaysis_vcf_compare
    - post_analaysis_phased_gene_variants
    - post_analaysis_lineage_occurance_table
    - post_analaysis_heatmap
    - post_analaysis_depth_histogram_figure
        - plots constructed depth data to varying plots
    - post_analaysis_dotplot
        - creates dotplot from mummer data
* remove_PCR_duplicates.smk
    - remove_PCR_duplicates
        - Uses Picard to remove PCR duplicates
    - bam_coverage
        - Uses ????? to determine mapping coverage at every position in genome. #TODO: add correct software
* run_trim.smk
    - run_trim
        - Trims reads of adapters and poor quality regions using BBduk
* samtools.smk
    - convert_sam2bam
        - Converts SAM alignment files to BAM files using Samtools
    - sort_bam
        - Sorts BAM files using Samtools
    - index_bam
        - Indexes BAM files using Samtools

## Options

## Usage

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/DGE-co-culture-workflow.git` or `git remote add -f upstream https://github.com/snakemake-workflows/DGE-co-culture-workflow.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.

### Step 7: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.
