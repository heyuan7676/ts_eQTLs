#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


### Collect results from all methods

ml R/3.5.1

for tau in 1000 500 100 50 20 10
do
	for seed in {1..10}
	do
		echo "Collecting results for tau${tau}_seed${seed}"
		Rscript simulation/compare_methods.R  -t ${tau} -s ${seed}
	done
done
