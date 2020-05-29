#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --no-requeue
#SBATCH -p quickdebug

source ./GLOBAL_VAR.sh

module load python/2.7

#outdir=${datadir}/input_pairs_fitModel
#mkdir -p ${outdir}


r2=1
#ith="$1"
#python filter_SNP_groupSNPs.py ${r2} ${suffix} ${ith}

python filter_SNPs_TagSNPs.py ${r2} ${suffix}

