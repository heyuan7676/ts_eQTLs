import os
import sys
import numpy as np
import pandas as pd

from readin_X import readin_X
from readin import readin
import pdb


r2 = '1' 
FDR = 0.05

fmDir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph'
ll_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL'
prefix       = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_%s' % str(r2)
LDprefix     = '_LD1'

FMfn = 'SparseMF_coph_%s_topPair_K30_a11_l110' % prefix.replace(r2, '0.2')
#FMfn  = 'SparseMF_coph_%s_topPair_K25_a125_l15000' % prefix
LMfn= '%s%s_Loadings_beta_BH_corrected_alpha%s' % (FMfn, LDprefix, str(FDR))
#LMfn = '%s%s_Loadings_projection' % (FMfn, LDprefix)

bg_cluster_id = 0


inputdatadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs'
inputdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel'
inputdatafn  = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance.txt'

pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets'


allSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs'
SNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/cbset_pairs'
datasetName   = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.SNP_loc'

activeSNPdir     = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active'
sigSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_features'
activeSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active_features'
active_proportion = 0.0

