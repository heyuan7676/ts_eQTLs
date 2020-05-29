import os
import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
import os
import sys
sys.setrecursionlimit(10000)
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tissue_specific_eQTL/DownStream')

import pandas as pd
import numpy as np


from scipy.io import loadmat
from statsmodels.stats import multitest 

from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from collections import Counter
from scipy.stats import percentileofscore

from scipy.stats import pearsonr
from scipy.stats import ranksums
from scipy.stats import ttest_ind
from sklearn.preprocessing import scale 


from GLOBAL_VAR import *



def filter_TF_TPM(tf_tmp_med, group, tfNames, TPM):
        passed_tfi = []
        
        ## filter out TFs with median TPM < 1 in matched tissues      
        # deal with names incosistency in TF files
        if 0:
            extra_tfi = list(set(tfNames) - set(tf_tmp_med.index))
            realNames_tfi   = [x.split('::') for x in extra_tfi]
            for i in range(len(realNames_tfi)):
                tfi = realNames_tfi[i]
                tfi = np.intersect1d(tfi, tf_tmp_med.index)
                if len(tfi) < len(realNames_tfi[i]):
                    continue
                else:
                    tfs_pass = np.sum(tf_tmp_med.loc[tfi][Comp_tissues[int(group)]]>TPM, axis=1) == len(Comp_tissues[int(group)])
                    if np.sum(tfs_pass) == len(tfi):
                        passed_tfi.append(extra_tfi[i])
                    
        passed_tfi = []
            
        # other consistent names
        tf_inexp_dataset = np.intersect1d(tfNames, tf_tmp_med.index)
        tf_tmp_med = tf_tmp_med.loc[tf_inexp_dataset]
        tpm_above_idx = np.where(np.sum(tf_tmp_med[Comp_tissues[int(group)]]>TPM, axis=1) > len(Comp_tissues[int(group)])/2.0)[0]
        tfs_pass_tpm = np.array(tf_tmp_med.index)[tpm_above_idx]
        
        allPassed_tf = np.unique(list(tfs_pass_tpm) + list(passed_tfi))
        #assert len(tfs_pass_tpm)+ len(passed_tfi) == len(allPassed_tf)

        return allPassed_tf



def readin_enrichment(LMfn, featureName, folds = 100):
    all_OR = pd.DataFrame()
    
    for theTis in ['Shared'] + list(tissues):
        try:
            fntemps = [x for x in os.listdir(enrich_dir) if ((LMfn in x) and (featureName in x) and (theTis in x) and ('NULL' not in x))]
	    np.random.seed(1)
            fntemps = np.random.choice(fntemps, folds, replace = False)
            dat_tfbs = pd.read_csv('%s/%s' % (enrich_dir, fntemps[0]), sep='\t', index_col = 0, usecols=[0,3,4,5,6])
	    dat_tfbs = dat_tfbs[~dat_tfbs.index.duplicated()]
            for fntemp in fntemps[1:]:
                dat_tfbs_tp = pd.read_csv('%s/%s' % (enrich_dir, fntemp), sep='\t', index_col = 0, usecols=[0,3,4,5,6])
		dat_tfbs_tp = dat_tfbs_tp[~dat_tfbs_tp.index.duplicated()]
                index_tp = np.unique(list(dat_tfbs.index) + list(dat_tfbs_tp.index))
                dat_tfbs = dat_tfbs.loc[index_tp]
                dat_tfbs_tp = dat_tfbs_tp.loc[index_tp]
                for col in dat_tfbs.columns:
                    dat_tfbs.loc[np.isnan(dat_tfbs[col]),col] = 0
                    dat_tfbs_tp.loc[np.isnan(dat_tfbs_tp[col]),col] = 0
                dat_tfbs = dat_tfbs + dat_tfbs_tp
        except:
            print('TFBS for %s does not exist ' % theTis)
            continue

        dat_tfbs['ts_TFBS']    = [np.ceil(x/len(fntemps)) for x in dat_tfbs['ts_TFBS']]
        dat_tfbs['ts_nonTFBS'] = [np.ceil(x/len(fntemps)) for x in dat_tfbs['ts_nonTFBS']]
        #dat_tfbs['background_TFBS']    = [np.ceil(x/len(fntemps)) for x in dat_tfbs['background_TFBS']]
        #dat_tfbs['background_nonTFBS'] = [np.ceil(x/len(fntemps)) for x in dat_tfbs['background_nonTFBS']]

    
        or_pv = dat_tfbs.apply(lambda x: fisher_exact([[x[2], x[3]], [x[0], x[1]]], alternative='greater'), axis=1)
        dat_tfbs['OR_background'] = [x[0] for x in or_pv]
        dat_tfbs['PV_background'] = [x[1] for x in or_pv]
        
        #left_tfs = filter_TF_TPM(tf_tmp_med, theTis, np.unique(dat_tfbs.index), TPM = TPM)
        #dat_tfbs = dat_tfbs.loc[left_tfs]
	dat_tfbs.index = [x.upper() for x in dat_tfbs.index]
	dat_tfbs['TF'] = dat_tfbs.index
	# restrict to single TFs
	single_tfs = np.intersect1d(dat_tfbs['TF'], tf_tmp_med.index)
	dat_tfbs = dat_tfbs.loc[single_tfs]
	if theTis == 'Shared':
		dat_tfbs['TF_TPM'] = dat_tfbs.apply(lambda x: np.mean(tf_tmp_med.loc[x['TF']]), axis=1)
	else:
		dat_tfbs['TF_TPM'] = dat_tfbs.apply(lambda x: np.mean(tf_tmp_med.loc[x['TF']][theTis]), axis=1)
        dat_tfbs['group'] = [theTis] * len(dat_tfbs)    

        all_OR = all_OR.append(dat_tfbs)
            
        print(theTis, dat_tfbs.shape)

    #all_OR = all_OR[all_OR['ts_TFBS'] + all_OR['background_TFBS'] > 10]
    #all_OR['pvalue2'] = all_OR.apply(lambda x: chi2_contingency([[x['ts_TFBS'], x['ts_nonTFBS']], [x['background_TFBS'], x['background_nonTFBS']]])[1], axis=1)
    
    ## GO score
    tfi_score = []       
    for i in range(len(all_OR)):
        rowi = all_OR.iloc[i]
        score_i = score_go(rowi['TF'], rowi['group'])
        #score_i = -1
        tfi_score.append(score_i)    
    all_OR['score'] = tfi_score

    all_OR['BH_pvalue'] = multitest.multipletests(all_OR['PV_background'], method = 'fdr_bh')[1]
    #all_OR['BH_pvalue_2'] = multitest.multipletests(all_OR['pvalue2'], method = 'fdr_bh')[1]
    all_OR = all_OR.sort_values('PV_background')
    dat_sig = all_OR[all_OR['BH_pvalue'] < 0.05]
    print('Using %s, for %s, Among %d(%d) tested TFs, %d(%d) are significantly enriched at FDR 0.05' % (LMfn_method, feature, len(all_OR), len(set(all_OR.index)), len(dat_sig), len(set(dat_sig['TF']))))

    return all_OR


def score_go(gi, ti):
    if ti == 'Shared':
	ts = tissues
    else:
    	ts = np.intersect1d(ti, geneset_weights.columns)
    if len(ts) == 0:
        return -1
    if gi not in gene_gos.keys():
        return -1
    try:
        #goi = gene_gos[gene_annotation[gi][0]]
        goi = gene_gos[gi]
        maxscore = np.nanmax(np.array(geneset_weights.loc[np.array(goi)][ts]).ravel())
    except:
        return 0
    return maxscore


def readin_gene_gos():
	### GO terms and the genes
	go_terms = {}
	gene_gos = {}
	GSEA_fn = open('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/c5.bp.v6.2.symbols.gmt.txt', 'r')
	for l in GSEA_fn.readlines():
		go = l.rstrip().split('\t')[0]
    		genes = ' '.join(l.rstrip().split('\t')[2:])
    		go_terms[go] = genes
    		for gi in l.rstrip().split('\t')[2:]:
        		try:
            			gene_gos[gi].append(go)
        		except:
            			gene_gos[gi] = [go]

	return gene_gos



def readin_go_weights():
        go_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/weighted_GO_terms/'
        fn = 'SetTestWeights_C5.BP_RNA.and.protein.txt'
        go_weights_df = pd.read_csv('%s/%s' % (go_dir, fn), sep='\t')
        
        go_weights_df = go_weights_df.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)), axis=0)

        go_weights_df.columns = ['Adipose_Subcutaneous', 'Adrenal_Gland', 'appendix', 'bone marrow',
                        'Breast_Mammary_Tissue', 'Brain_Cortex', 'Uterus', 'Colon_Sigmoid',
                        'duodenum', 'endometrium', 'epididymis', 'Esophagus_Mucosa', 'fallopian tube',
                        'gallbladder', 'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung',
                        'lymph node', 'ovary', 'Pancreas', 'parathyroid gland', 'placenta', 'prostate',
                        'rectum', 'salivary gland', 'seminal vesicle', 'Muscle_Skeletal', 'Skin_Not_Sun_Exposed_Suprapubic',
                        'Small_Intestine_Terminal_Ileum', 'smooth muscle', 'Spleen', 'Stomcah',
                        'Testis', 'Thyroid', 'tonsil', 'urinary bladder']
        go_weights_df['Esophagus_Muscularis'] = go_weights_df['Esophagus_Mucosa'].copy()
	go_weights_df['Whole_Blood'] = go_weights_df['Spleen'].copy()
                        
        go_weights_df = go_weights_df[np.intersect1d(go_weights_df.columns, tissues)]
        return go_weights_df


def readin_X_flashr():
    X = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/flashr/flashr_fm.txt', sep='\t')
    X = X / np.max(abs(X), axis=0)
    comp_tissues = []
    tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
    tissues = np.array(tissues[0])

    comp_tissues.append(tissues)
    for x in np.array(X.transpose())[1:, :]:
        comp_tissues.append(tissues[np.where(abs(x)  > 0.5)])
    return [X, comp_tissues, tissues]





if __name__ == '__main__':

	enrich_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features/0815_noRestriction'
	tf_tmp_med = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/GTEx_v8_median_gene_tpm.txt', sep='\t', index_col = 0)
	tf_tmp_med = tf_tmp_med.sort_values(tf_tmp_med.columns[0],  ascending = False)
	tf_tmp_med = tf_tmp_med.drop_duplicates()

	LMfn_method = 'Thresholding'
	folds = 2
 
	if LMfn_method == 'spMF':	
        	[X, Comp_tissues, tissues] = readin_X(FMfn)
	elif LMfn_method == 'flashr':
        	LMfn = 'flashr_Loadings_beta_BH_alpha0.05_corrected'
        	[X, Comp_tissues, tissues] = readin_X_flashr()
        elif LMfn_method == 'Thresholding':
                [X, Comp_tissues, tissues] = readin_X(FMfn)
                LMfn = 'Thresholding'

	gene_gos = readin_gene_gos()
	geneset_weights = readin_go_weights()

	for feature in ['ROADMAP_7_Enh', 'ROADMAP_1_TssA']:
		ROADMAP_df_tsboth_oneside =  readin_enrichment(LMfn, feature, folds=folds)
		ROADMAP_df_tsboth_oneside.to_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/plots/Fig5_TF_enrichment_%s_%s_iter%d_noRes.txt' % (LMfn_method, feature, folds), sep='\t')

