import pandas as pd
import numpy as np
import os
import sys

from collections import Counter
import pdb
from GLOBAL_VAR import *

def each_tis_specific(tis, folds_change_to_top = 100):
        ## read in the pairs in credible set
        print(tis)
        pairsets = []
        ciseQTL_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate'
        f = open('%s/%s_95set_pairs.txt' % (ciseQTL_dir, tis), 'r')
        for l in f.readlines():
            if gene_anno[l.rstrip().split('\t')[0]] == 'protein_coding':
                pairsets.append(' '.join((l.rstrip().split('\t')[0], l.rstrip().split('\t')[1])))
        f.close()

        ## pairs that is in the credible set for this tissue
        tis_cb_pairs = np.array(pairsets)
        print('    Number of pairs in the credible set: %d' % (len(tis_cb_pairs)))
        
        ## read in pvalues for all pairs in this tissue
        pv_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/cbset_ciseQTL_results'
        fn = '%s/%s.v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.txt' % ( pv_dir, tis)
        pv_df = pd.read_csv(fn, usecols=[0,1,6], sep = '\t', index_col = [0,1], header=None)
        pv_df.index = [' '.join(x) for x in pv_df.index]

        ## index the pairs in credible set
        potential_ts_pairs = pv_df.loc[tis_cb_pairs]

        ## get the pairs that have p-values close to top eQTLs
        potential_ts_pairs['gene'] = [x.split(' ')[0] for x in potential_ts_pairs.index]
        
        top_cbpairs_pv = potential_ts_pairs.groupby('gene')[6].min()
        potential_ts_pairs['top_pv'] = np.array(top_cbpairs_pv.loc[np.array(potential_ts_pairs['gene'])])
        tspairs_pv_close_to_top = potential_ts_pairs[np.abs(potential_ts_pairs[6] / potential_ts_pairs['top_pv']) < folds_change_to_top]
        print('    Number of pairs with pvalue close to the top eQTL: %d' % len(tspairs_pv_close_to_top))
    
        return np.array(tspairs_pv_close_to_top.index)



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


   

if __name__ == '__main__':
	FOLDS = 100
	PROP = 0.5
	PVALUE = 0.001
	N1 = 5

    	gene_anno = {}
    	gene_annoFn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt', 'r')
    	for line in gene_annoFn.readlines():
        	g = line.split('\t')[0]
        	gene_anno[g] = line.split('\t')[6]


	tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
	tissues = np.array(tissues[0])

	tissue_groups =  get_tis_groups()

	## step 1
	df_tissues = {}
	for tis in tissues:
    		df_tissues[tis] = each_tis_specific(tis, FOLDS)


	eQTL_in_tis_groups = {}
	noeffect_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/revision/noEffect_eQTLs'
	for tis_index in range(len(tissue_groups)):
		tis_g = tissue_groups[tis_index]
		## step 2
    		tis_snps = []
    		for t in tis_g:
        		tis_snps.append(df_tissues[t])
    		tis_snps = [a for b in tis_snps for a in b]
    		tis_snps_df = pd.DataFrame.from_dict(Counter(tis_snps), orient='index')
    		tis_snps_pairs = np.array(tis_snps_df[tis_snps_df[0] >= len(tis_g) * PROP ].index)
    
		## step 3
   		no_tis_N = (len(tissues)-len(tis_g)) - N1
    		noeffect_pairs = np.load('%s/pv_%s_tisN_%d.npz' % (noeffect_dir, str(PVALUE), no_tis_N))
    		noeffect_pairs = noeffect_pairs['idle_pairs']
 
    		eQTL_in_tis_groups = np.intersect1d(tis_snps_pairs, noeffect_pairs)

		LMfn = 'ts_closeToTop_FOLDS%d_PROP%s_PVALUE%s_N1%d' % (FOLDS, str(PROP), str(PVALUE), N1)
		outfn = '%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, tis_index)
		x = pd.DataFrame(eQTL_in_tis_groups)
		x.to_csv(outfn, sep='\t', index = False, header=False)



	### shared
		
