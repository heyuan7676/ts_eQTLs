#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p skylake

ml python/2.7


#python 6_GSEA_Genesets.py
for g in {-1..27}
do
	bash 6_GSEA.sh ${g} 5
done

python 6_GSEA_save.py 5
