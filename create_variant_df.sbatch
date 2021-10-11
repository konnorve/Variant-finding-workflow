#!/usr/bin/env bash
# 
#SBATCH --job-name=var_df
#
# Number of tasks/cores for job
# SBATCH -n 10
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 48:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="Running Variant Finding Workflow"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
# SBATCH --cpus-per-task=1
#
#SBATCH --mem 250000
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate base

python3 workflow/scripts/all_calls_df_creation.py