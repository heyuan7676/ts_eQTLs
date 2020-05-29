#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p shared

ml python/2.7

idx="$1"
fn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets/match_random/batches_random/Shared_ts_ciseQTL_closeToTop_random_matched_idx${idx}.txt
if [[ -f "${fn}" ]]
then
	rm ${fn}
	#exit	
fi
echo ${idx}
python 2_generate_random_pair_from_bg_thresholding.py Shared ${idx}
