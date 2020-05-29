#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p skylake

ml python/2.7


#python 6_Geneset_GSEA.py
for g in {0..22}
do
	bash 6_GSEA.sh ${g} 5
done
