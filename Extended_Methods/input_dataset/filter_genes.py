import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
import pandas as pd
import numpy as np
import os
#from matplotlib import pyplot as plt
from collections import Counter
import pdb


def filterGenes(allGenes):
    ### Filter out non-functional genes
    uniqGenes = np.unique(allGenes)

    gene_anno = {}
    gene_annoFn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt', 'r')
    for line in gene_annoFn.readlines():
        g = line.split('\t')[0]
        if g in uniqGenes:
            gene_anno[g] = line.split('\t')[6]

   
    allGenes_types = [gene_anno[g] for g in allGenes]
    allGenes_types = np.array(allGenes_types)
    ### allGenes_types: 69% protein coding, 9% linc RNA, 9% anti-sense, 5% processed_pseudogene, rest are other 32 gene types. 
    quailfy_pairs_idx = np.where(allGenes_types == 'protein_coding')[0]
    print('After filtering for genes that dont produce protein, there are %d pairs left from %d' % (len(quailfy_pairs_idx), len(allGenes)))

    #plt.figure()
    #pd.DataFrame.from_dict(Counter(allGenes_types), orient='index').plot.bar()
    #plt.title('Gene types for gene in all SNP-gene pairs')
    #pd.DataFrame.from_dict(Counter(gene_anno.values()), orient='index').plot.bar()
    #plt.title('Gene types for unique gene in all SNP-gene pairs')
    #plt.show()
    #plt.savefig('%s/gene_types.pdf' % outdir)
    #plt.close()


    return quailfy_pairs_idx




### read in all pairs
#datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/cbset_per_tissue_W0.9/aggregate'
#pairs = pd.read_csv('%s/v8_cbset_95_allPairs.txt' % datadir, sep='\t')

datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate'
pairs = pd.read_csv('%s/v8_cbset_95_allPairs.txt' % datadir, sep='\t')

good_idx = filterGenes(np.array(pairs['Gene']))
filtered_pairs = pairs.iloc[good_idx]

outfn = '%s/v8_cbset_95_allPairs_filteredGenes.txt' % datadir
filtered_pairs.to_csv(outfn, sep='\t', index=False)

allSNPs = set(filtered_pairs['SNP'])
snpfn = open('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/SNPIDs.txt', 'w')
for snp in allSNPs:
	snpfn.write('%s\n' % snp)






