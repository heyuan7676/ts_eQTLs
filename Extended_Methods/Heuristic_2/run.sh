#!/bin/bash
#SBATCH --time 6:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=2

ml python/2.7

#python 1_Remove_compensate_pairs.py
#python 1_Compute_nonzero_sets.py
#python 2_annotate_pairs.py

#for i in {0..22}
#do
#	for j in {0..4}
#	do
#		sbatch 2_generate_random_pair_from_bg_thresholding.sh ${i} ${j}
#	done
#done


python 5_Compute_active_regions_thr.py
bash 5_extract_active_SNP_features_thr.sh 7_Enh
bash 5_extract_active_SNP_features_thr.sh 1_TssA

for g in {-1..27}
do
        bash 5_TFBS_snps_genes_thr.sh ${g}
done


#bash 4_extract_SNP_features_thr.sh
#bash 4_extract_random_SNP_features_5folds_thr.sh

#python SNPset_active_enrichment_v2.py 0
#python SNPset_active_enrichment_v2.py 1



for g in {-1..22}
do
        for f in ROADMAP_1_TssA ROADMAP_7_Enh
	do
		for IDX in 0 1
		do
                                sbatch 5_TFBS_enrichment_multiple_sampling_separate_thr.sh ${g} ${f} ${IDX}
		done
        done
done

