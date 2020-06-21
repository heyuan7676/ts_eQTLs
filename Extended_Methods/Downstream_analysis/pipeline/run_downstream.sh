#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1

method="$1"
cd ../${method}
source ./GLOBAL_VAR.sh

### map eQTLs to factors
bash ${scripts_dir}/0_loadings/fitL.sh ${method}
ml python/2.7
python ${scripts_dir}/1_Remove_compensate_pairs.py
python ${scripts_dir}/1_Compute_nonzero_sets.py


### generate random SNP-gene pairs 
python ${scripts_dir}/2_annotate_pairs.py
python ${scripts_dir}/2_savePairN.py

for i in {0..22}
do
        for j in {0..4}
        do
                sbatch ${scripts_dir}/2_generate_random_pair_from_bg.sh ${i} ${j}
        done
done


### perform enrichment analysis of chromHMM status
bash ${scripts_dir}/4_extract_SNP_features.sh
python ${scripts_dir}/SNPset_active_enrichment.py 1


### perform enrichment analysis of TF motifs
python ${scripts_dir}/5_Compute_active_regions.py
bash ${scripts_dir}/5_extract_active_SNP_features.sh 7_Enh
bash ${scripts_dir}/5_extract_active_SNP_features.sh 1_TssA
for g in {0..22}
do
        python ${scripts_dir}/5_TFBS_snps_genes.py ${g} ROADMAP_7_Enh
        python ${scripts_dir}/5_TFBS_snps_genes.py ${g} ROADMAP_1_TssA
done



