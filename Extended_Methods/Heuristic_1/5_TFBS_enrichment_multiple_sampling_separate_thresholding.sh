#!/bin/bash
#SBATCH --time 30:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p lrgmem

 
group="$1"
featureName=ROADMAP_7_Enh
ml python/2.7
ml parallel

source ./GLOBAL_VAR.sh

idx=0
fn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features/0815_noRestriction/TF_enrichment_ROADMAP_7_Enh_THR400_${group}_TFBS_Thresholding_torandom_idx${idx}.txt
if [ ! -f ${fn} ]; then
	echo $fn
	python 5_TFBS_enrichment_multiple_sampling_separate_noRes_thresholding.py ${group} ${featureName} ${idx}
fi


idx=1
fn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features/0815_noRestriction/TF_enrichment_ROADMAP_7_Enh_THR400_${group}_TFBS_Thresholding_torandom_idx${idx}.txt
if [ ! -f ${fn} ]; then
        echo $fn
        python 5_TFBS_enrichment_multiple_sampling_separate_noRes_thresholding.py ${group} ${featureName} ${idx}
fi
