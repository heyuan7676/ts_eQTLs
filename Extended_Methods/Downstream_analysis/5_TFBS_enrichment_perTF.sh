#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p shared

source ./GLOBAL_VAR.sh

g="$1"
idx="$2"
startidx=$((idx*60))
endidx=$(( $startidx + 59 ))
echo $startidx, $endidx

module load python/2.7
for (( idx=$startidx; idx<=$endidx; idx++ ))
do
	for f in ROADMAP_7_Enh ROADMAP_1_TssA
	do
		python ${scripts_dir}/5_TFBS_enrichment_perTF.py "$g" "$f" 0 "$idx"
		python ${scripts_dir}/5_TFBS_enrichment_perTF.py "$g" "$f" 1 "$idx"
	done
done

