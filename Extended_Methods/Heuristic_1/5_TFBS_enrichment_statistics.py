import os
import sys
import numpy as np
import pandas as pd

from readin_X import readin_X
import pdb

from GLOBAL_VAR import *
from scipy.stats import fisher_exact
from scipy import stats
from scipy.stats import ranksums

def readin_outliers(theGROUP):
    ## outliers
    outLier = dict()
    outLier_df = pd.read_csv('%s/%s_outlierPairs%sgroup%s.txt' % (pairdir, LMfn, randomPattern, theGROUP), sep='\t', header=None)
    outLier_df.columns = ['Gene', 'SNP']
    outLier[theGROUP] = outLier_df

    if background == 'shared':
    	outLier_df = pd.read_csv('%s/%s_outlierPairs%sgroup%s.txt' % (pairdir, LMfn, randomPattern, 0), sep='\t', header=None)
    	outLier_df.columns = ['Gene', 'SNP']
    else:
        outLier_df = pd.read_csv('%s/%s_outlierPairs_%sgroup%s.txt' % (pairdir, LMfn, 'random_matched_', theGROUP), sep='\t', header=None)
        outLier_df.columns = ['Gene', 'SNP']
    outLier[-1] = outLier_df.drop_duplicates()

    return outLier



def small_format(listoflist):
        list_pairs = [a for b in listoflist for a in b]
        pairs_df_new = pd.DataFrame([p.split(':') for p in list_pairs])
        pairs_df = pairs_df_new.copy()
        pairs_df.columns = ['Gene', 'SNP']
        pairs_df = pairs_df.set_index('SNP', drop = False)
	return pairs_df





def intersection_in_enhancer(ref_pair):
    ## eQTL in enhancer 
    outlier_fc = pd.read_csv('%s/%s_outlierSNPs%s%s_group%d.txt' % (sigSNPfeaturedir, LMfn, randomPattern, featureName.split('_')[0], theGROUP), sep='\t', index_col = 0)
    outlier_fc.columns = [x.replace('chr1_','') for x in outlier_fc.columns]
    existing_tissue = np.intersect1d(Comp_tissues[theGROUP], outlier_fc.columns)

    ref_feature = outlier_fc.loc[ref_pair['SNP']][existing_tissue]
    eQTL_in_feature = np.zeros(len(ref_feature))
    for tis in existing_tissue:
        eQTL_in_feature = eQTL_in_feature + np.array(['_'.join(featureName.split('_')[1:]) in xx for xx in np.array(ref_feature[tis])]) * 1
    ref_pair = ref_pair.iloc[np.where(eQTL_in_feature)[0]]

    if background == 'shared':
        outlier_fc = pd.read_csv('%s/%s_outlierSNPs%s%s_group%d.txt' % (sigSNPfeaturedir, LMfn, randomPattern, featureName.split('_')[0], 0), sep='\t', index_col = 0)
    else:
        outlier_fc = pd.read_csv('%s/%s_outlierSNPs_random_matched_%s_group%d.txt' % (sigSNPfeaturedir, LMfn, featureName.split('_')[0], theGROUP), sep='\t', index_col = 0)
    outlier_fc.columns = [x.replace('chr1_','') for x in outlier_fc.columns]

    return ref_pair





if __name__ == '__main__':
    featureName = 'ROADMAP_7_Enh'
    background = 'random'
    THR = 400
    randomPattern = "_"

    [X, Comp_tissues, tissues] = readin_X(FMfn)

    ## genes with SNPs in TFBS
    TF_allSNP_median = []
    TF_allGenes_median = []
    dist_median = []
    dist_enh_median = []

    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features'
    for theGROUP in range(23):
    	try:
		TF_allsnp_genes = np.load('%s/TFBS_allsnpgenes/TFBS_allsnpgenes_%s_group%d_THR%d_%s.npz' % (outdir, featureName, theGROUP, THR, FMfn))
    	except:
		print('TFBS for Loading %d does not exist' % theGROUP)
		continue

    	TF_allSNP = TF_allsnp_genes['TFBS_snps'][()]
    	TF_allGenes = TF_allsnp_genes['TFBS_genes'][()]

    	## outliers
    	outLier = readin_outliers(theGROUP)
    	ts_pair = outLier[theGROUP]
    	ts_pair_enh = intersection_in_enhancer(ts_pair)

	TF_allSNP_median.append([theGROUP, np.median([len(x) for x in TF_allSNP.values()])])
	TF_allGenes_median.append([theGROUP, np.median([len(x) for x in TF_allGenes.values()])])
	dist_median.append([theGROUP, np.median(ts_pair.groupby('Gene').size())])
	dist_enh_median.append([theGROUP, np.median(ts_pair_enh.groupby('Gene').size())])


    pd.DataFrame(TF_allSNP_median).to_csv('%s/TF_allSNP_median.txt' % outdir, sep='\t')
    pd.DataFrame(TF_allGenes_median).to_csv('%s/TF_allGenes_median.txt' % outdir, sep='\t')
    pd.DataFrame(dist_median).to_csv('%s/dist_median.txt' % outdir, sep='\t')
    pd.DataFrame(dist_enh_median).to_csv('%s/dist_enh_median.txt' % outdir, sep='\t')





