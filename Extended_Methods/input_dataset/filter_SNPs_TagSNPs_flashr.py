import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')

import pandas as pd
import numpy as np
import os

from scipy.stats.mstats import gmean
import pdb
import time





def readin_features(feature, pairs):
    Df = pd.read_csv('%s/%s.%s.txt' % (inputdir, data_prefix, feature), sep = '\t', index_col = [0,1], header=None)
    pairs_in_tuple = zip(pairs['Gene'], pairs['SNP'])
    index_df = Df.loc[pairs_in_tuple]
    index_df.index.names = ['Gene', 'SNP']
    #df = index_df.sort_index()
    return index_df





def index_topSNP_pergene(pairDf):
 
    min_pv_per_pair = pd.DataFrame(pairDf[pairDf.columns[2:]].apply(lambda x: np.min(x), axis=1))
    min_pv_per_pair['Gene'] = pairDf['Gene']

    min_pair_idx = min_pv_per_pair.groupby('Gene')[0].apply(lambda x: np.argmin(x))
    pair_min_pv = pairDf.iloc[min_pair_idx]

    return pair_min_pv[['Gene', 'SNP']]



if __name__ == '__main__':
    ### read in all pairs
    r2 = str(sys.argv[1])
    data_prefix = sys.argv[2]
    datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets'
    inputdir = '%s/input_pairs' % datadir
    outdir = '%s/input_pairs_fitModel' % datadir

    # save the data for pairs with top SNP per gene
    pairDf = pd.read_csv('%s/%s_filteredSNPs.LDblocks_%s_pvalue.txt' % (outdir, data_prefix, r2), sep='\t')
    top_pairs = index_topSNP_pergene(pairDf)

    slope_Df = readin_features('slope', top_pairs)
    se_Df    = readin_features('se', top_pairs)

    save_prefix = '%s_filteredSNPs.LDblocks_%s' % (data_prefix, r2)
    slope_Df.to_csv('%s/%s_topPair_slope_flashr.txt' % (outdir, save_prefix), sep = '\t')
    se_Df.to_csv('%s/%s_topPair_se_flashr.txt' % (outdir, save_prefix), sep = '\t')

