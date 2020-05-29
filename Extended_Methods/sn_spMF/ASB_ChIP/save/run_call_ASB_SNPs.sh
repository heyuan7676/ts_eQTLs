#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


for TF in HNF4A CTCF FOXA1 FOXA2 ATF3 JUND MAX ZBTB33
do
	for g in 5 8 10 12 15 20
	do
		sbatch call_ASB_SNPs.sh ${TF} $g
	done
done
