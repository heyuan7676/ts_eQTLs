#!/bin/bash
#SBATCH --time 32:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p shared

#K="$1"
#alpha1="$2"
#lambda1="$3"


K=30
alpha1=49
lambda1=490
ml R/3.5.1
Rscript run_MF.R -k $K -a ${alpha1} -l ${lambda1} -t 50
