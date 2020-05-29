#!/bin/bash

FDR=0.05
pairdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets

prefix=v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete
LDprefix=_LD1
datasetName=${prefix}.SNP_loc

inputdatadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs

pairdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets

activeSNPdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active
sigSNPfeaturedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_features/
activeSNPfeaturedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active_features

allSNPfeaturedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs/
#SNPfeaturedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/cbset_pairs/

#FMfn=SparseMF_coph_${prefix}_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110
#FMfn=SparseMF_coph_${prefix}_filteredSNPs.LDblocks_1_topPair_K25_a125_l15000
#LMfn=${FMfn}${LDprefix}_Loadings_beta_BH_corrected_alpha${FDR}
#LMfn=${FMfn}${LDprefix}_Loadings_projection

FMfn=Thresholding_Refined
LMfn=ts_closeToTop_FOLDS100_PROP0.5_PVALUE0.001_N15
