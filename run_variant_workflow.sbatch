#!/usr/bin/env bash
#SBATCH --job-name=varsmk
#SBATCH --time 20-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o slurm.out
#SBATCH -e slurm.err

source activate snakemake

snakemake --snakefile workflow/Snakefile --unlock

snakemake \
    --snakefile workflow/Snakefile \
    --use-conda \
    --rerun-incomplete \
    --cluster "sbatch -p {resources.partition} -t {resources.time} -n {resources.ntasks} --mem={resources.mem} -o {resources.output} -e {resources.error}" \
    --jobs 10 \
    --keep-going