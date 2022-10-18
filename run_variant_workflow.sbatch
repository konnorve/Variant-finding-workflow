#!/usr/bin/env bash
#SBATCH --job-name=varsmk
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate snakemake

snakemake --snakefile workflow/Snakefile --unlock # --lint

# run job:
snakemake \
    --snakefile workflow/Snakefile \
    --use-conda \
    --conda-frontend conda \
    --latency-wait 60 \
    --cluster "sbatch -p newnodes,sched_mit_hill,sched_any -t 12:00 -N 1 -n 16 --mem=64G -o logs/smk_slurm/%j_slurm_output.txt -e logs/smk_slurm/%j_slurm_error.txt" \
    --jobs 10
    
# get report once job is finished
# snakemake --snakefile $SNAKEFILE --report report.html