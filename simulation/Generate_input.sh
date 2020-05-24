#!/bin/bash
#SBATCH --time 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p skylake



seed="$1"

for tau in 1000 100 10 1 20 500 50 5
do

	ml R/3.5.1
	Rscript simulation/Generate_input.R -t ${tau} -s ${seed}
done
