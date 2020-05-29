
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from scipy.io import loadmat

import os
import sys
sys.setrecursionlimit(10000)
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')


import matplotlib 
from matplotlib import pyplot as plt
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'

from matplotlib import cm
import seaborn as sns
sns.set_style('whitegrid')

from collections import Counter
from scipy.stats import ttest_ind
from itertools import combinations
from scipy.stats import spearmanr
from statsmodels.stats import multitest
from scipy.stats import spearmanr

from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial.distance import pdist

from sklearn.cluster import KMeans
from scipy.stats import probplot


# In[2]:


if 1:
    gene_anno = {}
    gene_annoFn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt', 'r')
    for line in gene_annoFn.readlines():
        g = line.split('\t')[0]
        gene_anno[g] = line.split('\t')[6]

        
tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
tissues = np.array(tissues[0])


# In[4]:


allpairs_pv_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs'
fn = '%s/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.pvalue.txt' % allpairs_pv_dir

allpairs_pv = pd.read_csv(fn, sep='\t', index_col= [0,1], header = None)
allpairs = np.array([' '.join(x) for x in allpairs_pv.index])
allpairs_pv.index = allpairs
allpairs_pv.columns = tissues


# In[3]:


def derive_ts_eQTL_one_tissue(noeffect_in_other_tis_pairs, tis, folds_change_to_top = 100):
    
        ## read in the pairs in credible set
        pairsets = []
        ciseQTL_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate'
        f = open('%s/%s_95set_pairs.txt' % (ciseQTL_dir, tis), 'r')
        for l in f.readlines():
            if gene_anno[l.rstrip().split('\t')[0]] == 'protein_coding':
                pairsets.append(' '.join((l.rstrip().split('\t')[0], l.rstrip().split('\t')[1])))
        f.close()

        ## pairs that is in the credible set for this tissue
        tis_cb_pairs = np.intersect1d(np.array(pairsets), noeffect_in_other_tis_pairs)
        print('    Number of pairs with high pvalues in other %d tissues, and in the credible set: %d' % (non_sig_tissue_N, len(tis_cb_pairs)))
        
        ## read in pvalues for all pairs in this tissue
        pv_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/cbset_ciseQTL_results'
        fn = '%s/%s.v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.txt' % ( pv_dir, tis)
        pv_df = pd.read_csv(fn, usecols=[0,1,6], sep = '\t', index_col = [0,1], header=None)
        pv_df.index = [' '.join(x) for x in pv_df.index]

        ## index the pairs in credible set
        cbpairs_pv = pv_df.loc[np.array(pairsets)]
        potential_ts_pairs = pv_df.loc[tis_cb_pairs]

        ## get the pairs that have p-values close to top eQTLs
        cbpairs_pv['gene'] = [x.split(' ')[0] for x in cbpairs_pv.index]       
        potential_ts_pairs['gene'] = [x.split(' ')[0] for x in potential_ts_pairs.index]
        
        top_cbpairs_pv = cbpairs_pv.groupby('gene')[6].min()
        potential_ts_pairs['top_pv'] = np.array(top_cbpairs_pv.loc[np.array(potential_ts_pairs['gene'])])
        tspairs_pv_close_to_top = potential_ts_pairs[np.abs(potential_ts_pairs[6] / potential_ts_pairs['top_pv']) < folds_change_to_top]
        print('    Number of pairs with pvalue close to the top eQTL: %d' % len(tspairs_pv_close_to_top))
    
        ## save
        outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/ciseQTL_cut'
        fn = open('%s/%s_ts_ciseQTL_closeToTop.txt' % (outdir,tis), 'w')
        for eQTL in np.array(tspairs_pv_close_to_top.index):
            [snpi, genei]= eQTL.split(' ')
            fn.write('%s\t%s\n' % (snpi, genei))
        fn.close()

        ## save at pairdir for downstream analysis
        pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets'
        fn = open('%s/%s_ts_ciseQTL_closeToTop.txt' % (pairdir,tis), 'w')
        for eQTL in np.array(tspairs_pv_close_to_top.index):
            [snpi, genei]= eQTL.split(' ')
            fn.write('%s\t%s\n' % (snpi, genei))
        fn.close()
        
        ## save the pvalues for these pairs -- for drawing examples
        ts_eQTL_pv_df = allpairs_pv.loc[np.array(tspairs_pv_close_to_top.index)]
        ts_eQTL_pv_df['top_eQTL_pv'] = tspairs_pv_close_to_top['top_pv']
        ts_eQTL_pv_df.to_csv('%s/%s_ts_ciseQTL_closeToTop_pvalues.txt' % (outdir,tis), sep='\t')

        return [len(tis_cb_pairs), len(ts_eQTL_pv_df)]



# In[4]:


non_sig_pv_thr = 0.001
non_sig_tissue_N = len(tissues) - 5
folds_change_to_top = 100

def derive_ts_eQTL():

    noeffect_in_other_tis = np.where(np.sum(allpairs_pv > non_sig_pv_thr, axis=1) > non_sig_tissue_N)[0]
    noeffect_in_other_tis_pairs = allpairs[noeffect_in_other_tis]
    
    filter1, filter2 = [], []
    for t in tissues:
        print(t)
        [A, B] = derive_ts_eQTL_one_tissue(noeffect_in_other_tis_pairs, t)
        filter1.append(A)
        filter2.append(B)
        print('')
        
    return [filter1, filter2]
        


# In[39]:


[N1,N2] = derive_ts_eQTL()


# In[46]:


df_plot = pd.DataFrame({"N1": N1, "N2": N2, "tissues": tissues})
df_plot = df_plot.melt(id_vars='tissues')


plt.figure(figsize = (12,6))
sns.barplot(x = 'tissues', y = 'value', hue='variable', data = df_plot)
plt.xticks(rotation = 90)
plt.show()
plt.close()


# In[49]:


savedir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/plots/'

df_plot = pd.DataFrame({"N1": N1, "N2": N2, "tissues": tissues})
df_plot.to_csv('%s/Fig2_ciseQTL_naive_closeToTop.txt' % savedir, sep='\t')


# In[ ]:


### shared

noeffect_in_other_tis = np.where(np.sum(allpairs_pv > non_sig_pv_thr, axis=1) > non_sig_tissue_N)[0]


# In[50]:


allpairs_pv.head()


# In[20]:


### sanity check

tis = 'Whole_Blood'
outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/ciseQTL_cut'
tp = pd.read_csv('%s/%s_ts_ciseQTL_closeToTop_pvalues.txt' % (outdir,tis), sep='\t', index_col = 0)


# In[21]:


## read in pvalues for all pairs in this tissue
pv_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/cbset_ciseQTL_results'
fn = '%s/%s.v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.txt' % ( pv_dir, tis)
pv_df = pd.read_csv(fn, usecols=[0,1,6], sep = '\t', index_col = [0,1], header=None)
pv_df.columns = ['pv']
pv_df.index = [' '.join(x) for x in pv_df.index]
pv_df['Gene'] = [x.split(' ')[0] for x in pv_df.index]
pv_df['SNP_loc'] = [int(x.split(' ')[1].split('_')[1]) for x in pv_df.index]


# In[47]:


sns.set(font_scale=2)
sns.set_style('whitegrid')

for i in np.random.choice(range(len(tp)), 10, replace = False)[:1]:
    plt.figure(figsize=(20,5))
    pv_vec = tp.iloc[i]
    
    ## p-values for this pair across 49 tissues
    plt.subplot(121)
    plt.scatter(range(49), -np.log10(np.array(pv_vec)[:49]))
    plt.scatter(list(tissues).index(tis), -np.log10(pv_vec[tis]), color = 'red')
    plt.xticks(range(49), tissues, rotation = 90, size = 8)
    plt.title('%s \n across 49 tissues' % tp.index[i])
    plt.ylabel('-log10(p-value)')

    [gi, snpi] = tp.index[i].split(' ')
    snpi_loc = int(snpi.split('_')[1])
    tp_i = pv_df[pv_df['Gene'] == gi]
    starti = np.min(tp_i['SNP_loc'])
    x_loc = np.array(tp_i['SNP_loc']) - starti
    logpv = -np.log10(np.array(tp_i['pv']))
    
    ## p-values for SNPs for this gene in this tissue
    plt.subplot(122)
    plt.scatter(x_loc, logpv)
    plt.xticks([], [])
    plt.scatter(snpi_loc - starti, -np.log10(pv_vec[tis]), color = 'red')
    plt.xlabel('Chromosome location')
    plt.title('all variants for %s\n in %s' % (gi, tis))
    plt.ylabel('-log10(p-value)')
    plt.show()
    plt.close()


# In[13]:


counts = Counter(allpairs)
counts_values = np.array(list(counts.values()))
counts_eQTL = np.array(list(counts.keys()))


# In[5]:


ciseQTL_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate'

def readin_pairSet():
    pairsets= {}
    for k in tissues:
        print(k)
        pairsets[k] = []
        #f = open('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn.replace('_uniq',''), k), 'r')
        f = open('%s/%s_95set_pairs.txt' % (ciseQTL_dir, k), 'r')
        for l in f.readlines():
            if gene_anno[l.rstrip().split('\t')[0]] == 'protein_coding':
                pairsets[k].append(' '.join((l.rstrip().split('\t')[0], l.rstrip().split('\t')[1])))
        f.close()

    return pairsets


tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
tissues = np.array(tissues[0])
pairsets = readin_pairSet()

allpairs = [a for b in pairsets.values() for a in b]
counts = Counter(allpairs)
counts_values = np.array(list(counts.values()))
counts_eQTL = np.array(list(counts.keys()))


# In[15]:


higher = 5
shared_eQTL = counts_eQTL[np.where(counts_values > higher)[0]]


outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/ciseQTL_cut'

## shared
fn = open('%s/Shared_ts_ciseQTL_closeToTop.txt' % outdir, 'w')
for eQTL in shared_eQTL:
    [genei, snpi]= eQTL.split(' ')
    fn.write('%s\t%s\n' % (genei, snpi))
fn.close()



pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets'

## shared
fn = open('%s/Shared_ts_ciseQTL_closeToTop.txt' % pairdir, 'w')
for eQTL in shared_eQTL:
    [genei, snpi]= eQTL.split(' ')
    fn.write('%s\t%s\n' % (genei, snpi))
fn.close()


# In[18]:


shared_eQTL

