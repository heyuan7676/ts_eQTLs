import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')

import pandas as pd
import numpy as np
import os

from scipy.stats.mstats import gmean
from scipy.io import savemat
import pdb
import time



def readinPvalueDf():
    ### obtain pvalue and pairs information
    if 0:
        Df = pd.read_csv('%s/%s.pvalue.txt' % (inputdir, data_prefix), sep = '\t', index_col = [0,1], header=None, nrows = None)
        gmean_pvalue = Df.apply(lambda x: gmean(x[~np.isnan(x)]), axis=1)
        gmean_pvalue.to_csv('%s/%s.pvalue.geman.txt' % (inputdir, data_prefix), sep='\t')
    else:
        gmean_pvalue = pd.read_csv('%s/%s.pvalue.geman.txt' % (inputdir, data_prefix), sep='\t', index_col = [0, 1], header=None)
        gmean_pvalue.columns = range(gmean_pvalue.shape[1])
        Pairs = np.array(gmean_pvalue.index)
        genes = [p[0] for p in Pairs]
        SNPs  = [p[1] for p in Pairs]
        gmean_pvalue.index = pd.MultiIndex.from_arrays([genes, SNPs])

    return gmean_pvalue


def getTagSNPs(gene_snp_pvalue):

    genes = [x[0] for x in gene_snp_pvalue.index]
    snps  = [x[1] for x in gene_snp_pvalue.index]
    pairs = [' '.join((x[0], x[1])) for x in gene_snp_pvalue.index]
    gene_snp_pairs = pd.DataFrame({'pairs':pairs, 'snps':snps, 'genes': genes})
    df = gene_snp_pairs.set_index('pairs')

    ## read in
    leftout_pairs = list()
    for ith in range(39):
	print ith
    	save_prefix = '%s_filteredSNPs.LDblocks_%s_ith_%d' % (data_prefix, r2, ith)
    	outFn   = '%s/npz_ith/%s.npz' % (outdir, save_prefix)
    	ld_pairs = np.load(outFn)
    	ld_pairs = ld_pairs['arr_0'].item()
    	leftout_pairs.append([a for b in list(ld_pairs.values()) for a in b])

    leftout_pairs = [a for b in leftout_pairs for a in b]
    tag_df = df.loc[np.array(list(set(df.index) - set(leftout_pairs)))]
    tag_pairs = np.sort(np.array(tag_df.index))
    print 'Use tag SNPs to represent each LD blocks: %d pairs left' % (len(tag_pairs))

    return tag_pairs





def readin_features(feature, pairs):
    Df = pd.read_csv('%s/%s.%s.txt' % (inputdir, data_prefix, feature), sep = '\t', index_col = [0,1], header=None)
    pairs_in_tuple = [(p.split(' ')[0], p.split(' ')[1]) for p in pairs]
    index_df = Df.loc[pairs_in_tuple]
    index_df.index.names = ['Gene', 'SNP']
    #df = index_df.sort_index()
    return index_df





def index_oneSNP_pergene(pairDf):
    pairDf.set_index('Gene', drop = False, inplace=True)
    np.random.seed(0)
    idx = np.random.permutation(range(len(pairDf)))
    pairDf = pairDf.iloc[idx]
    top_pairs = pairDf.groupby(pairDf.index).first()
    top_pairs = np.array(top_pairs.apply(lambda x: ' '.join((x['Gene'], x['SNP'])), axis=1))

    top_pair_fn = open('%s/%s_topPair.txt' % (outdir, save_prefix), 'w')
    for pr in top_pairs:
        top_pair_fn.write('%s\t%s\n' % (pr.split(' ')[0], pr.split(' ')[1]))
    top_pair_fn.close()

    return top_pairs



if __name__ == '__main__':
    ### read in all pairs
    r2 = str(sys.argv[1])
    data_prefix = sys.argv[2]
    datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets'
    inputdir = '%s/input_pairs' % datadir
    outdir = '%s/input_pairs_fitModel' % datadir

    s = time.time()
    gene_snp_pvalue =  readinPvalueDf()
    print 'Reading in the geman palues takes %f min' % ((time.time() - s) / 60.0)

    s = time.time()
    tag_pairs = getTagSNPs(gene_snp_pvalue)
    print 'Merging results of tag SNPs takes %f min' % ((time.time() - s) / 60.0)

    ## save
    save_prefix = '%s_filteredSNPs.LDblocks_%s' % (data_prefix, r2)
    outFn   = '%s/%s' % (outdir, save_prefix)
    fn = open('%s.txt' % outFn, 'w')
    for p in tag_pairs:
        fn.write('%s\n' % p)
    fn.close()

    # save the data for pairs after filtering SNPs
    slope_Df = readin_features('slope', tag_pairs)
    slope_Df.to_csv('%s/%s_slope.txt' % (outdir, save_prefix), sep = '\t')

    seDf     = readin_features('se', tag_pairs)
    seDf.to_csv('%s/%s_se.txt' % (outdir, save_prefix), sep = '\t')

    pvalueDf = readin_features('pvalue', tag_pairs)
    pvalueDf.to_csv('%s/%s_pvalue.txt' % (outdir, save_prefix), sep = '\t')


    # save the data for pairs with one SNP per gene
    if 1:
    	pairDf = pd.read_csv('%s.txt' % outFn, sep=' ', header=None)
    	pairDf.columns = ['Gene', 'SNP']

    	top_pairs = index_oneSNP_pergene(pairDf)
    	slope_Df = readin_features('slope', top_pairs)
    	se_Df    = readin_features('se', top_pairs)
    	slope_Df.to_csv('%s/%s_topPair_slope.txt' % (outdir, save_prefix), sep = '\t')
    	se_Df.to_csv('%s/%s_topPair_se.txt' % (outdir, save_prefix), sep = '\t')

    	# save in mat
    	result = dict()
    	result['X'] = np.array(slope_Df)
    	result['SD'] = np.array(se_Df)
   	savemat('%s/%s_topPair.mat' % (outdir, save_prefix), result)


