#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p debug

method="$1"
g="$2"
idx="$3"

cd ../${method}
scripts_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/revision/Other_methods/pipeline/scripts

ml python/2.7
echo ${method}


startidx=$((idx*100))
endidx=$(( $startidx + 99 ))
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


