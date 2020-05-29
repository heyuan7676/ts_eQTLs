import os
import sys
import pandas as pd
import numpy as np
from GLOBAL_VAR import *
from collections import Counter



def append_df_change_index(gi_dat, mm, nn, mm_l, nn_l):
    try:
    	df_temp = gi_dat.loc[mm_l, nn_l]
    except:
	return gi_dat

    index_array = [[int(mm)] * len(df_temp), [int(nn)] * len(df_temp)]
    df_temp.index = index_array
    df = gi_dat.append(df_temp)
    df = df.sort_index()
    return df


def match_for_randomPairs(tis, randomBatch, allGenes):
    pair_Fn = '%s/match_random/%s_ts_ciseQTL_closeToTop.txt' % (pairdir, tis)
    pairs_df    = pd.read_csv(pair_Fn, sep='\t', low_memory=False)
    pairs_df.columns = ['Genes', 'SNPs', 'dis', 'dis_cat', 'maf', 'maf_cat', 'pair']

    gtex_allpairs_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/allPairs/maf_dis/split_by_gene'
    output_Fn = open('%s/match_random/batches_random/%s_ts_ciseQTL_closeToTop_random_matched_idx%d.txt' % (pairdir, tis, randomBatch), 'w')

    not_matched_genes = []
    gi = ""

    for pair_slice in pairs_df.groupby('Genes'):
        attemps = 0
	print(pair_slice[0])
        snps_for_this_gene = np.array(pair_slice[1][['dis_cat', 'maf_cat']])
	allGenes = allGenes[allGenes!=gi]
	    
        while True:
            attemps += 1
            if attemps > 1000:
                not_matched_genes.append(pair_slice[0])
                print('%s: not matchable' % pair_slice[0])
                break

            gi = np.random.choice(allGenes, 1)[0]
            if gi == pair_slice[0]:
                continue

    
            fn = '%s/Gene_%s_bg.txt' % (gtex_allpairs_dir, gi)
	    gi_dat = pd.read_csv(fn, sep=' ', header=None)
	    gi_dat = gi_dat.set_index([4,5]).sort_index()
           
	    true_sets = Counter([':'.join((str(int(x[0])), str(int(x[1])))) for x in snps_for_this_gene])
	    pool_sets = Counter([':'.join((str(int(x[0])),str(int(x[1])))) for x in gi_dat.index])

	    ### make sure there are enough SNP-genes pairs for matching 
	    gi_dat.columns = ['gene', 'snp', 'dis', 'maf']
	    random_snp_vec = []
	    nomatched = 0
	    for token in true_sets.keys():
		[mm, nn] = token.split(':')
		if true_sets[token] > pool_sets[token]:
			try:
				# relax MAF match criteria
				nn_l = max(int(nn) - 1, 0)
			except:
				nomatched = 1
				#print(token)
				break
			nn_h = min(int(nn) + 1, 10)
			mm_l = max(int(mm) - 1, 0)
			mm_h = min(int(mm) + 1, 50)
			gi_dat = append_df_change_index(gi_dat, mm, nn, mm_l, nn_l)
			gi_dat = append_df_change_index(gi_dat, mm, nn, mm_l, nn_h)
			gi_dat = append_df_change_index(gi_dat, mm, nn, mm_h, nn_l)
			gi_dat = append_df_change_index(gi_dat, mm, nn, mm_h, nn_h)

			pool_sets_updated = Counter([':'.join((str(x[0]),str(x[1]))) for x in gi_dat.index])
			if token not in pool_sets_updated.keys():
				nomatched = 1
			elif pool_sets_updated[token] < true_sets[token]:
				nomatched = 1

			if nomatched:
				#print('%s has no properly matched pairs' % gi)
				break

		snp_pool_updated = gi_dat.loc[int(mm), int(nn)]
		snps_i = np.random.choice(np.array(snp_pool_updated['snp']),true_sets[token], replace = False)
		random_snp_vec.append(snps_i)

	    if nomatched:
		continue
		
	    #print('%s got the matched pairs' % gi)	
	    random_snp_vec = [a for b in random_snp_vec for a in b]
	    assert len(random_snp_vec) == len(snps_for_this_gene)
	
	    for snpi in random_snp_vec:
		output_Fn.write('%s\t%s\n' % (gi, snpi))

	    break
    
    output_Fn.close()



if __name__ == '__main__':
	allGenes = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/allGenes.txt', sep='\t', header=None)
	allGenes = np.array(allGenes[0])

	tis = sys.argv[1] ## should be tissues
	idx = 1 * int(sys.argv[2]) #0 --> 0 to 49
	for k in range(1):
                #try:
		#	output_Fn = pd.read_csv('%s/match_random/batches_random/%s_ts_ciseQTL_closeToTop_random_matched_idx%d.txt' % (pairdir, tis, idx), sep='\t', nrows = 10, header =None)
                #except:
		match_for_randomPairs(tis, idx, allGenes)
		print('Idx%d finished for %s' % (idx, tis))
		idx += 1



