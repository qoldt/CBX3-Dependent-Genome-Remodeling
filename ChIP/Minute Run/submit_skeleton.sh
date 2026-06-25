#!/bin/bash

#SBATCH --job-name=minute
#SBATCH --output=logs/%x-%j.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --partition=long
#SBATCH --time=72:00:00

export TMPDIR=/data/cephfs-1/home/users/newmana_c/scratch/tmp
export LOGDIR=logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
mkdir -p $LOGDIR


eval "$($(which conda) shell.bash hook)"
conda activate minuteman

set -x

minute run --scheduler greedy --keep-going --retries 3 --jobs 60 --rerun-incomplete --cluster-config /data/cephfs-1/home/users/newmana_c/work/Minute/cluster.yaml --cluster 'sbatch --time {cluster.time} --partition {cluster.partition} -c {cluster.cpus} --ntasks-per-node=1 -e logs_slurm/{cluster.jobname}.err -o logs_slurm/{cluster.jobname}.out -J {cluster.jobname}'
