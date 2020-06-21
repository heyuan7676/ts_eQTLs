#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -p skylake

ml gcc/5.5.0
ml R/3.5.1


idx="$1"
FMfn="$2"
Rscript lm.R ${idx} ${FMfn}




