#!/bin/bash

scripts_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/Revision_Zscore/scripts

allSNPfeaturedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs/
inputdatadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs

route_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/Revision_Zscore
ll_dir=${route_dir}/LL
pairdir=${route_dir}/downstream/pairSets
activeSNPdir=${route_dir}/downstream/SNP/SNPset_active
sigSNPfeaturedir=${route_dir}/downstream/SNP/SNPset_features/
activeSNPfeaturedir=${route_dir}/downstream/SNP/SNPset_active_features

enrichment_dir=${route_dir}/downstream/enrichmentTest/

gsea_fn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/c5.bp.v6.2.symbols.gmt.txt
gsea_dir=${enrichment_dir}/GSEA

FMfn=PMD_cv1
LMfn=${FMfn}_Loadings_beta_BH_alpha0.05_corrected
