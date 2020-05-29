

## Generate ts-eQTLs and shared-eQTLs
python ts_ciseQTL_naive_revision.py
python shared_ciseQTL_naive_revision.py

## annotate the eQTLs
python 2_annotate_pairs.py


## generate matched random pairs
for i in {-1..27}
do
	for j in {0..4}
	do
		sbatch 2_generate_random_pair_from_bg_thresholding.sh ${i} ${j}
	done
done


## extract DNase & ROADMAP features for the eQTLs and random pairs
bash 4_extract_random_SNP_features_5folds_thr.sh
bash 4_extract_SNP_features_thr.sh


## compute enrichment of eSNPs / random SNPs in DNase & ROADMAP features
python SNPset_active_enrichment_v2_thresholding.py 0 
python SNPset_active_enrichment_v2_thresholding.py 1


## extract open regions for the FMfn
python 5_Compute_active_regions_thr.py


## annotate the TF motif for eQTLs in Enh / TssA
bash 5_extract_active_SNP_features_thr.sh

## extract eGenes with TF motif & compute enrichment in TF motif
for i in {-1..27}
do
	for j in ROADMAP_7_Enh ROADMAP_1_TssA
	do	
		sbatch 5_TFBS_snps_genes_thr.sh ${i} ${j}
	done
done


## compute enrichment in TF motif
for i in {-1..27}
do
        for j in ROADMAP_7_Enh ROADMAP_1_TssA
        do
                sbatch 5_TFBS_snps_genes_thr.sh ${i} ${j}
        done
done





