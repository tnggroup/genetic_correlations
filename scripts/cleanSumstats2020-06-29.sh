#!/bin/bash -l

#SBATCH -n 8

#SBATCH -t 1:00:00

/mnt/lustre/groups/ukbiobank/sumstats/scripts/cleanSumstats.py $1 $2
