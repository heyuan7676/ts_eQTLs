import os
import sys
import numpy as np
import pandas as pd

from readin_X import readin_X
import pdb

from GLOBAL_VAR import *


#### collect SNPs in TFBS
#### 1) get the TFs where there are at least 1 SNP in the TFBS
#### to exclude those super-conservative TFs where SNPs never happen on TFBS

#### 2) in next step, restrict to genes with at least 1 SNP in the TFBS

def readin_geneannotation():
    ### read in gene annotation
    gene_annotation = {}
    gene_annotation_fn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt','r')
    for l in gene_annotation_fn.readlines():
	l = l.rstrip().split('\t')
        ensemID = l[0].replace('"','')
        gene_annotation[ensemID] = [l[1], l[2]]
    return gene_annotation





if __name__ == '__main__':
    THR = int(sys.argv[1])
    group = int(sys.argv[2])
    featureName = sys.argv[3]  # ROADMAP_7_Enh

    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features/TFBS_allsnpgenes'
    [X, Comp_tissues, tissues] = readin_X(FMfn)

    gene_annotation = readin_geneannotation()
    allGenes = pd.read_csv('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/allGenes.txt', sep='\t',header=None)
    allGenes.columns = ['Gene']
    allGenes['chr'] = [gene_annotation[g][0] for g in allGenes['Gene']]
    allGenes['tss'] = [gene_annotation[g][1] for g in allGenes['Gene']]
    allGenes['tss_m1mb'] = [np.max([int(x)-1000000, 0]) for x in allGenes['tss']]
    allGenes['tss_p1mb'] = [int(x)+1000000 for x in allGenes['tss']]


    if 1:
	#### collect TFBS in all active regions
	TFBS_snps, TFBS_genes = {}, {}

	print featureName, group
	#SNP_TFBS = pd.read_csv('%s/%s_Active_SNPset_%s_TFmotif_group%d.txt' % (activeSNPfeaturedir, FMfn, featureName, group), sep = '\t', index_col = 0)	
	SNP_TFBS = pd.read_csv('%s/%s_Active_SNPset_%s_TFmotif_Esophagus_Mucosa.txt' % (activeSNPfeaturedir, FMfn, featureName), sep = '\t', index_col = 0) 
	if len(SNP_TFBS) == 0:
		sys.exit()
	snps      = np.array(SNP_TFBS.index)

	for tfi in SNP_TFBS.columns:
		TFBS_snps[tfi.replace('chr1_','')]  = np.unique(snps[np.where(SNP_TFBS[tfi] > THR)[0]])
		snps_df = pd.DataFrame({"SNP": TFBS_snps[tfi.replace('chr1_','')] })
		snps_df['SNP_chr'] = [x.split('_')[0] for x in snps_df['SNP']]
		snps_df['SNP_loc'] = [int(x.split('_')[1]) for x in snps_df['SNP']]

		TFBS_genes_tfi = []
		for snps_chr in set(snps_df['SNP_chr']):
			gene_chr = allGenes[allGenes['chr'] == snps_chr]
			snp_loc_matrix = np.reshape(np.repeat(np.array(snps_df['SNP_loc']), len(gene_chr)), [len(snps_df), len(gene_chr)])
			gene_tss_m1mb_matrix = np.reshape(np.repeat(np.array(gene_chr['tss_m1mb']), len(snps_df)), [len(gene_chr), len(snps_df)]).transpose()
			gene_tss_p1mb_matrix = np.reshape(np.repeat(np.array(gene_chr['tss_p1mb']), len(snps_df)), [len(gene_chr), len(snps_df)]).transpose()

			qualify_pair_idx = np.where(((gene_tss_m1mb_matrix - snp_loc_matrix) < 0) * ((gene_tss_p1mb_matrix - snp_loc_matrix) > 0))
			snps_inTSS  = np.array(snps_df['SNP'])[qualify_pair_idx[0]]
			genes_inTSS = np.array(gene_chr['Gene'])[qualify_pair_idx[1]]
			TFBS_genes_tfi.append(np.unique(genes_inTSS))
			
		tfi = tfi.replace('chr1_','')
		TFBS_genes[tfi] = np.unique([a for b in TFBS_genes_tfi for a in b])
		print('For %s, there are %d genes and %d snps' % (tfi, len(TFBS_snps[tfi]), len(TFBS_genes[tfi])))

        outfn = '%s/TFBS_allsnpgenes_%s_group%d_THR%d_%s_Esophagus_Mucosa' % (outdir, featureName, group, THR, FMfn)
	#outfn = '%s/TFBS_allsnpgenes_%s_group%d_THR%d_%s' % (outdir, featureName, group, THR, FMfn)
        print(outfn)
        np.savez(outfn, TFBS_snps = TFBS_snps, TFBS_genes = TFBS_genes)







