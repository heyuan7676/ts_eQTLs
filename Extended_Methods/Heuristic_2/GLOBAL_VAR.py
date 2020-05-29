import os
import sys
import numpy as np
import pandas as pd

import pdb


r2 = '1' 
FDR = 0.05

fmDir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph'
ll_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL'
prefix       = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_%s' % str(r2)
LDprefix     = '_LD1'

#FMfn = 'SparseMF_coph_%s_topPair_K30_a11_l110' % prefix.replace(r2, '0.2')
#FMfn  = 'SparseMF_coph_%s_topPair_K25_a125_l15000' % prefix
#LMfn= '%s%s_Loadings_beta_BH_corrected_alpha%s' % (FMfn, LDprefix, str(FDR))
#LMfn = '%s%s_Loadings_projection' % (FMfn, LDprefix)

FMfn = 'Thresholding_Refined'


if 1:
        FOLDS = 100
        PROP = 0.5
        PVALUE = 0.001
        N1 = 5
LMfn = 'ts_closeToTop_FOLDS%d_PROP%s_PVALUE%s_N1%d' % (FOLDS, str(PROP), str(PVALUE), N1)

bg_cluster_id = 0


inputdatadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs'
inputdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel'
inputdatafn  = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance.txt'

pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets'
#pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets_0907'


allSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs'
SNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/cbset_pairs'
datasetName   = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.SNP_loc'

activeSNPdir     = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active'
sigSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_features'
activeSNPfeaturedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/SNP/SNPset_active_features'
active_proportion = 0.0

gsea_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/GSEA/'
gsea_file_used = 'c5.bp.v6.2.symbols.gmt.txt'

def get_tis_groups():
        tissue_groups = [[x for x in tissues if 'Adipose' in x],
                 ['Adrenal_Gland'],
                [x for x in tissues if 'Artery' in x],
                [x for x in tissues if 'Brain' in x],
                 ['Cells_EBV-transformed_lymphocytes'],
                 ['Cells_Cultured_fibroblasts'],
                [x for x in tissues if 'Colon' in x],
                [x for x in tissues if 'Esophagus' in x],
                [x for x in tissues if 'Heart' in x],
                 ['Kidney_Cortex'],
                 ['Liver'],
                 ['Lung'],
                 ['Minor_Salivary_Gland'],
                 ['Muscle_Skeletal'],
                 ['Nerve_Tibial'],
                 ['Ovary'],
                 ['Pancreas'],
                 ['Pituitary'],
                 ['Prostate'],
                [x for x in tissues if 'Skin' in x],
                ['Small_Intestine_Terminal_Ileum'],
                ['Spleen'],
                ['Stomach'],
                ['Testis'],
                ['Thyroid'],
                ['Uterus'],
                ['Vagina'],
                ['Whole_Blood']]
        return tissue_groups


tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
tissues = np.array(tissues[0])

Comp_tissues = get_tis_groups()
