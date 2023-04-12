#!/usr/bin/env bash
#SBATCH --job-name=var-smk
#SBATCH --time 5-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o var-smk.out
#SBATCH -e var-smk.err

source activate snakemake

snakemake --profile profile
