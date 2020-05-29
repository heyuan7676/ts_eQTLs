import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')

import pandas as pd
import numpy as np
import os

from scipy.stats.mstats import gmean
import time


def readinPvalueDf():
    ### obtain pvalue and pairs information
    Df = pd.read_csv('%s/%s.pvalue.txt' % (inputdir, data_prefix), sep = '\t', index_col = [0,1], header=None, nrows = 1000)
    gmean_pvalue = Df.apply(lambda x: gmean(x[~np.isnan(x)]), axis=1)
    gmean_pvalue.to_csv('%s/%s.pvalue.geman.txt' % (inputdir, data_prefix), sep='\t')
    return gmean_pvalue




if __name__ == '__main__':
    ### read in all pairs
    r2 = str(sys.argv[1])
    data_prefix = sys.argv[2]
    datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets'
    inputdir = '%s/input_pairs' % datadir

    s = time.time()
    gene_snp_pvalue =  readinPvalueDf()
    print 'Reading in the geman palues takes %f min' % ((time.time() - s) / 60.0)




