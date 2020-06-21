#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p skylake


method="$1"
cd ../${method}
source ./GLOBAL_VAR.sh

ml python/2.7


python ${scripts_dir}/6_GSEA_Genesets.py
for g in {0..22}
do
	bash ${scripts_dir}/6_GSEA.sh ${g} 5
done

python ${scripts_dir}/6_GSEA_save.py 6
