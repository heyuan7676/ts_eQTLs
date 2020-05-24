#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


### Collect results from all methods

tau="$1"
seed="$2"


ml R/3.5.1
Rscript simulation/compare_methods.R  -t ${tau} -s ${seed}
