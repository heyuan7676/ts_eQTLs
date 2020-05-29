#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -p debug

idsDir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate
rm -f ${idsDir}/v8_cbset_95_allPairs.txt
rm -f ${idsDir}/temp
rm -f ${idsDir}/temp_sorted

for f in `ls ${idsDir}/*_95set_pairs.txt`
do
	awk '{print $1,$2}' ${f} | sed 's/ /	/g' >> ${idsDir}/temp
done

cat ${idsDir}/temp | sort  > ${idsDir}/temp_sorted

echo "Gene	SNP" > ${idsDir}/v8_cbset_95_allPairs.txt
cat ${idsDir}/temp_sorted | uniq >> ${idsDir}/v8_cbset_95_allPairs.txt

module load python/2.7
python filter_genes.py

