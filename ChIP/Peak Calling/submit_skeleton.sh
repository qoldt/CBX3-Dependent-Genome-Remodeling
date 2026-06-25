#!/bin/bash 
#SBATCH --job-name=macs3
#SBATCH --output=logs/%x-%j.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --partition=medium
#SBATCH --time=24:00:00
 
export TMPDIR=/data/cephfs-1/home/users/newmana_c/scratch/tmp
export LOGDIR=logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
mkdir -p $LOGDIR 

eval "$($(which conda) shell.bash hook)"
conda activate minuteman
set -x
snakemake --profile=cubi-v1 -j 48 -k -p --restart-times=2 
