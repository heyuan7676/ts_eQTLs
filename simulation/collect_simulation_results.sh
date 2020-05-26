#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p skylake


tau="$1"
seed="$2"

ml R/3.5.1
Rscript simulation/tune_parameters.R -O simulation/output/tau${tau}_seed${seed}/
