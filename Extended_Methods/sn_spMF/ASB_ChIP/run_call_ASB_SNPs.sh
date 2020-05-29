#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


for TF in HNF4A CTCF
do
	for g in 6 8 10 12 15
	do
		sbatch run_ASB.sh ${TF} $g 1
	done
done
