#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --no-requeue
#SBATCH -p quickdebug

source ./GLOBAL_VAR.sh

module load python/2.7

r2=1

python filter_SNPs_TagSNPs_flashr.py ${r2} ${suffix}

