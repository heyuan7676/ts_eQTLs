#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

for f in `ls /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/*bed | xargs -n1 basename`; do echo $f;  sbatch annotate_TF_chromHMM.sh $f; done

#for f in `ls /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/*bed | xargs -n1 basename `; do f=${f%.bed}; sbatch extract_SNP_features.sh ${f} ""; sbatch extract_SNP_features.sh ${f} "_random_matched"; done

