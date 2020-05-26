import os
import sys
import numpy as np
import pandas as pd
from scipy.io import loadmat

import pdb


FDR = 0.05
prefix       = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1'
LDprefix     = '_LD1'

inputdatadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs'
inputdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel'
inputdatafn  = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance.txt'

allSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs'
datasetName   = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.SNP_loc'

gtex_allpairs_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/allPairs/maf_dis/split_by_gene'


route_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/Revision'

ll_dir = '%s/LL' % route_dir
pairdir = '%s/downstream/pairSets' % route_dir
activeSNPdir     = '%s/downstream/SNP/SNPset_active' % route_dir
sigSNPfeaturedir = '%s/downstream/SNP/SNPset_features' % route_dir
activeSNPfeaturedir = '%s/downstream/SNP/SNPset_active_features' % route_dir
active_proportion = 0.0

enrichment_dir = '%s/downstream/enrichmentTest/' % route_dir
chromHMM_dir = '%s/active_regions' % enrichment_dir
gsea_fn = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/c5.bp.v6.2.symbols.gmt.txt'
gsea_dir = '%s/GSEA' % enrichment_dir

plot_save_dir = '%s/plots/' % route_dir

tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
tissues = np.array(tissues[0])


## read in Factor matrix
def get_factor_tissues(FMfn, tissues):
    X = pd.read_csv('%s.txt' % FMfn, sep='\t')
    X = np.array(X / np.max(abs(X), axis=0))
    X = X.transpose()

    Comp_tissues = []
    for x in np.array(X):
        Comp_tissues.append(tissues[np.where(abs(x)  > 0.01)])
    return [X, Comp_tissues]
