from GLOBAL_VAR import *
from collections import OrderedDict



def readin_pairSet():
    snpsets, pairsets= {}, {}
    #for k in range(K):
    for k in tissues:
	pairsets[k] = []
	snpsets[k]  = []
        #f = open('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn.replace('_uniq',''), k), 'r')
	#f = open('%s/%s_95set_pairs.txt' % (ciseQTL_dir, k), 'r')
	#f = open('%s/flashr_Loadings_beta_BH_alpha0.05_corrected_outlierGenes_group%d.txt' % (pairdir, k), 'r')
	f = open('%s/%s_ts_ciseQTL_closeToTop.txt' % (pairdir, k), 'r')
	for l in f.readlines():
		pairsets[k].append(' '.join((l.rstrip().split('\t')[0], l.rstrip().split('\t')[1])))
		snpsets[k].append(l.rstrip().split('\t')[1])
        f.close()

    return [pairsets, snpsets]



def readin_TSS():
    if 0:
        dis_tss = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance_notused.txt', sep='\t', header=None)
	genes = dis_tss[0]
	snps  = dis_tss[1]
        del dis_tss[0]
        del dis_tss[1]
	dis = np.abs(dis_tss.apply(lambda x: np.unique(x)[0], axis=1))
	dis_df = pd.DataFrame({"Genes":genes, "SNPs":snps, "dis": dis})
	dis_df.to_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance.txt', sep='\t', index = False)
    else:
	dis_df = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.tss_distance.txt', sep='\t')

    genes = dis_df['Genes']
    snps  = dis_df['SNPs']
    dis_df.index = [' '.join(p) for p in zip(genes,snps)]

    return dis_df



def readin_SNP_MAF():
    SNP_AF  = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/v8_cbset_95_SNPs_AF.txt', sep='\t', header=None, index_col = 0)
    SNP_AF['maf'] = SNP_AF[1].apply(lambda x: np.min([x, 1-x]))
    SNP_AF = SNP_AF[~SNP_AF.index.duplicated(keep='first')]

    return SNP_AF




if __name__ == '__main__':

    ciseQTL_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/caviar_output_GTEx_LD/aggregate'
    [X, Comp_tissues, tissues] = readin_X(FMfn)
    #tissues = list(tissues) + ['Shared']
    tissues = ['Shared']
    K = len(Comp_tissues)

    # two metrices to annotate
    SNP_AF = readin_SNP_MAF()
    SNP_AF['maf_cat'] = [int(x*100)/5 for x in SNP_AF['maf']]

    # 
    DIS_TSS = readin_TSS()
    DIS_TSS['dis_cat'] = [int(x/20000) for x in DIS_TSS['dis']]

    ### get real SNP sets
    [Pair_sets, SNP_sets] = readin_pairSet()
    #for k in range(K):
    for k in tissues:
	print(k)
	outlier_SNP_AF = SNP_AF.loc[np.unique(list(SNP_sets[k]))]
	outlier_SNP_AF = outlier_SNP_AF[['maf', 'maf_cat']]
	outlier_pair_DIS = DIS_TSS.loc[np.unique(list(Pair_sets[k]))]

    	df_pairs = outlier_pair_DIS.merge(outlier_SNP_AF, left_on = 'SNPs', right_index = True)
    	df_pairs['pair'] = df_pairs.index
	df_pairs.to_csv('%s/match_random/%s_ts_ciseQTL_closeToTop.txt' % (pairdir, k), sep = '\t', index = False)
	#df_pairs.to_csv('%s/match_random/flashr_Loadings_beta_BH_alpha0.05_corrected_outlierGenes_group%d_annotated.txt' % (pairdir,k), sep='\t', index=False)
	#df_pairs.to_csv('%s/match_random/%s_outlierPairs_group%d_annotated.txt' % (pairdir, LMfn, k), sep='\t', index=False)
	#df_pairs.to_csv('%s/match_random/ciseQTL_%s.txt' % (pairdir, k), sep='\t', index=False)

