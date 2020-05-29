#!/bin/bash
#SBATCH --time 6:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p unlimited

tfi_idx="$1" ## 567 - 572
feature=ROADMAP_1_TssA

ml python/2.7


source ./GLOBAL_VAR.sh

for idx in 10 11
do
	python 5_TFBS_enrichment_multiple_sampling_separate_noRes_oneTF_thresholding.py Shared ${feature} ${idx} ${tfi_idx}
done

