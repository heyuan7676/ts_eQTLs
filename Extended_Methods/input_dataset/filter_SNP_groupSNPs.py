import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')

import pandas as pd
import numpy as np
import os

import networkx as nx
import pickle
from scipy.stats.mstats import gmean
import pdb
import time

def deriveLDblocks(r2_thr):
	LD_block_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes'
	fn = '%s/plink_r2_%s.ld.graphs.pickle' % (LD_block_dir, r2_thr)
	#fn = '%s/temp.graphs.pickle' % (LD_block_dir)
	with open(fn, 'rb') as input:
		graph = pickle.load(input)
	comp = list()
	for c in nx.connected_components(graph):
		comp.append(list(c))
	print 'For r2>%s, There are %d LD blocks' % (str(r2_thr), len(comp))
	return comp


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



def filterSNPs(gene_snp_pvalue, blocks, ith = None):

    genes = [x[0] for x in gene_snp_pvalue.index]
    snps  = [x[1] for x in gene_snp_pvalue.index]
    pairs = [' '.join((x[0], x[1])) for x in gene_snp_pvalue.index]
    gene_snp_pairs = pd.DataFrame({'pairs':pairs, 'snps':snps, 'genes': genes})
    df = gene_snp_pairs.set_index('snps')

    ### obatin the LD blocks
    ld_pairs = dict()

    i = 0
    gap = 10000
    if ith is not None:
	startidx = ith * gap
	endidx   = (ith + 1) * gap
	if startidx > len(blocks):
	    sys.exit('ith exceeds the number of blocks')
	if endidx > len(blocks):
	    endidx = len(blocks) 
	blocks = blocks[startidx:endidx]
	print 'Blocks - %d to %d' % (startidx, endidx)

    for c in blocks:
        pool = np.array(c)
	pool_in_df = np.intersect1d(pool, df.index)
	if len(pool_in_df) <= 1:
		continue
	else:
            df_for_this_pool = df.loc[pool_in_df]
        for df_slice in df_for_this_pool.groupby('genes'):
	    if len(df_slice[1]) == 1:
		continue
            snps = df_slice[1].index
            tp =  snps[np.argmin([gene_snp_pvalue.loc[df_slice[0]].loc[snpii][0] for snpii in snps])]
            snps = list(snps)
            snps.remove(tp)
            leftout_pairs = np.array(df_slice[1].loc[np.array(snps)]['pairs'])
            ld_pairs[df_slice[1].loc[tp]['pairs']] = leftout_pairs

    ## save
    save_prefix = '%s_filteredSNPs.LDblocks_%s_ith_%d' % (data_prefix, r2, ith)
    outFn   = '%s/%s' % (outdir, save_prefix)
    np.savez(outFn, ld_pairs)



if __name__ == '__main__':
    ### read in all pairs
    r2 = str(sys.argv[1])
    data_prefix = sys.argv[2]
    ith = int(sys.argv[3])
    datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets'
    inputdir = '%s/input_pairs' % datadir
    outdir = '%s/input_pairs_fitModel/npz_ith' % datadir

    s = time.time()
    graph = deriveLDblocks(r2)
    print 'Reading in the graphs takes %f min' % ((time.time() - s) / 60.0)

    s = time.time()
    gene_snp_pvalue =  readinPvalueDf()
    print 'Reading in the geman palues takes %f min' % ((time.time() - s) / 60.0)

    s = time.time()
    filterSNPs(gene_snp_pvalue, graph, ith = ith)
    print 'Filtering SNPs takes %f min' % ((time.time() - s) / 60.0)

