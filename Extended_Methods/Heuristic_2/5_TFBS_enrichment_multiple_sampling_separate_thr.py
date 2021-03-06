import os
import sys
import numpy as np
import pandas as pd

import pdb

from GLOBAL_VAR import *
from scipy.stats import fisher_exact
from scipy import stats
from scipy.stats import ranksums
import time




def readin_outliers(theGROUP):
    ## outliers
    outLier = dict()
    outLier_df = pd.read_csv('%s/%s_outlierPairs%sgroup%s.txt' % (pairdir, LMfn, randomPattern, theGROUP), sep=' ', header=None)
    outLier_df.columns = ['Gene', 'SNP']
    outLier[theGROUP] = outLier_df

    outLier_df = pd.DataFrame()
    for i in [idx, idx+3]:
        outLier_dfi = pd.read_csv('%s/match_random/batches_random/%s_outlierPairs_%sgroup%s_idx%d.txt' % (pairdir, LMfn, 'random_matched_', theGROUP, i), sep='\t', header=None)
        outLier_dfi.columns = ['Gene', 'SNP']
        outLier_df = outLier_df.append(outLier_dfi)

    outLier[-100] = outLier_df

    return outLier



def small_format(listoflist):
        list_pairs = [a for b in listoflist for a in b]
        pairs_df_new = pd.DataFrame([p.split(':') for p in list_pairs])
        pairs_df = pairs_df_new.copy()
        pairs_df.columns = ['Gene', 'SNP', 'GeneID']
        pairs_df = pairs_df.set_index('SNP', drop = False)
	return pairs_df



def match_SNPset_size_Gene(ref_pairs, other_pairs):
        ## deal with duplicated genes in bg

        dist1 = ref_pairs.groupby('Gene').size().sort_values(ascending = True)
        dist2 = other_pairs.groupby('Gene').size().sort_values(ascending = True)

        dist2_permuted = pd.Series()
        for k in np.sort(np.unique(dist2))[::-1]:
                x = dist2[dist2 == k]
                x = x.iloc[np.random.choice(range(len(x)), len(x), replace = False)]
                dist2_permuted = dist2_permuted.append(x)
        dist2 = dist2_permuted.copy()

        if len(dist1) > len(dist2):
                dist1 = dist1[:len(dist2)]

        pool1 = ref_pairs.groupby('Gene')['SNP'].apply(list)
        pool2 = other_pairs.groupby('Gene')['SNP'].apply(list)

        ref_genes = ref_pairs.set_index('Gene', drop = False)
        other_genes = other_pairs.set_index('Gene', drop = False)

        matched_pairs_ref, matched_pairs_other = [], []

        for rowi  in range(len(dist1)):
                gi_ref = dist1.index[rowi]
                gi = np.argmin(np.abs(dist2 - dist1.iloc[rowi]))
                #gi = dist2.index[rowi]
		#print(dist1.iloc[rowi] - dist2.loc[gi])

                if dist1.iloc[rowi] < dist2.loc[gi]:
                        tp = pool1[gi_ref]
                        matched_pairs_ref.append([':'.join((np.unique(ref_genes.loc[gi_ref]['Gene'])[0], s, str(gi_ref))) for s in tp])
                        tp = np.random.choice(pool2[gi], dist1.iloc[rowi], replace = False)
                        matched_pairs_other.append([':'.join((np.unique(other_genes.loc[gi]['Gene'])[0], s, str(gi))) for s in tp])
                else:
                        tp = np.random.choice(pool1[gi_ref], dist2.loc[gi], replace = False)
                        matched_pairs_ref.append([':'.join((np.unique(ref_genes.loc[gi_ref]['Gene'])[0], s, str(gi_ref))) for s in tp])
                        tp = pool2[gi]
                        matched_pairs_other.append([':'.join((np.unique(other_genes.loc[gi]['Gene'])[0], s, str(gi))) for s in tp])
                dist2 = dist2.drop(gi)


        other_pairs_new = small_format(matched_pairs_other)
        ref_pairs_new   = small_format(matched_pairs_ref)

        return [ref_pairs_new, other_pairs_new]








def separate_pairs(ts_pair, background_pair, countAs = 'NONE'):
        # ts-pairs 
        if countAs == 'NONE':
	    ts_genes = np.array(list(set(ts_pair['Gene']) - set(background_pair['Gene'])))
            ts_pair = ts_pair.set_index('Gene', drop = False).loc[ts_genes].set_index('SNP', drop=False)

	    background_genes = np.array(list(set(background_pair['Gene']) - set(ts_pair['Gene'])))
            background_pair = background_pair.set_index('Gene', drop = False).loc[background_genes].set_index('SNP', drop=False)
            
        elif countAs == 'ts':
	    background_genes = np.array(list(set(background_pair['Gene']) - set(ts_pair['Gene'])))
            background_pair = background_pair.set_index('Gene', drop = False).loc[background_genes].set_index('SNP', drop=False)
            
	# deal with duplicated genes in bg
        xs = list(background_pair['Gene'])
        unique_gene_bg = [i for i, n in enumerate(xs) if i==0 or n != xs[i-1]]
        gene_add_one = np.zeros(len(background_pair))
        gene_add_one[np.array(unique_gene_bg)] = 1
        unique_gene_bg = np.cumsum(gene_add_one)
        background_pair['GeneID'] = unique_gene_bg

	xs = list(ts_pair['Gene'])
        unique_gene = [i for i, n in enumerate(xs) if i==0 or n != xs[i-1]]
        gene_add_one = np.zeros(len(ts_pair))
        gene_add_one[np.array(unique_gene)] = 1
        unique_gene = np.cumsum(gene_add_one)
        ts_pair['GeneID'] = unique_gene

        return [ts_pair, background_pair]








def find_asscio_pairs_for_snps_in_TFBS_compare_rr(ts_pair, background_pair, TF_genes, TF_snps):

        gene1 = np.intersect1d(ts_pair['Gene'], TF_genes)
        gene2 = np.intersect1d(background_pair['Gene'], TF_genes)

        ts_pair = ts_pair.set_index('Gene', drop=False).loc[gene1]
        background_pair = background_pair.set_index('Gene', drop = False).loc[gene2]

        [ts_pair, background_pair]     = match_SNPset_size_Gene(ts_pair, background_pair)

        snp_in_tfbs = np.intersect1d(ts_pair['SNP'], TF_snps)
        snp_in_tfbs_bg = np.intersect1d(background_pair['SNP'], TF_snps)

        ts_pair['In'] = False
        background_pair['In'] = False

        ts_pair = ts_pair.set_index('SNP', drop = False)
        ts_pair.loc[snp_in_tfbs, 'In'] = True

        background_pair = background_pair.set_index('SNP', drop = False)
        background_pair.loc[snp_in_tfbs_bg, 'In'] = True


        ts_gene_prop = ts_pair.groupby('Gene').apply(lambda x: np.mean(x['In']))
        bg_gene_prop = background_pair.groupby('Gene').apply(lambda x: np.mean(x['In']))

        A = np.sum(ts_gene_prop > 0)
        B = len(ts_gene_prop) - A
        C = np.sum(bg_gene_prop > 0)
        D = len(bg_gene_prop) - C

        return [A, B,C,D]




def contig_test(a, b, c, d):
        [a,b,c,d] = [float(a), float(b), float(c), float(d)]
        if c == 0:
                c += 1
                d += 1
                a += 1
                b += 1
	rs = fisher_exact([[a,b], [c,d]], alternative = 'greater')
        return [(a/b)/(c/d), rs[1]]




def collect_results(pair_set1, pair_set2):
    ## collect results
    OR1, PV1, OR2, PV2, TF = [], [], [], [], []
    N1, N2, N3, N4, N5, N6 = [], [], [], [], [], []
    for tfi in np.sort(list(TF_allSNP.keys())):
        if len(set(TF_allSNP[tfi])) == 0:
                continue
        [A, B, C, D] = find_asscio_pairs_for_snps_in_TFBS_compare_rr(pair_set1, pair_set2, TF_allGenes[tfi],TF_allSNP[tfi])

        tfi = tfi.split('_')[1]
        if (B==0) or (D == 0):
            continue
        if contig_test(A, B, C, D)[1] < 1:
            print(theTis, tfi, A, B, C,D, contig_test(A, B, C, D))

        [oddsratio, pv] = contig_test(A, B, C, D)
        OR1.append(oddsratio)
        PV1.append(pv)

        TF.append(tfi)
        N1.append(A)
        N2.append(B)
        N3.append(C)
        N4.append(D)

    result = pd.DataFrame({"TF": TF,  "OR_background": OR1, "PV_background": PV1, "ts_TFBS": N1, "ts_nonTFBS": N2, "background_TFBS": N3, "background_nonTFBS": N4})
    result = result.set_index('TF').sort_values('PV_background')

    return result






if __name__ == '__main__':
    theTis = int(sys.argv[1])
    featureName = sys.argv[2]
    idx = int(sys.argv[3])

    THR = 400
    randomPattern = "_"

    print("%s - idx%d" % (theTis,idx ))

    ## genes with SNPs in TFBS
    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features'
    try:
	TF_allsnp_genes = np.load('%s/TFBS_allsnpgenes/TFBS_allsnpgenes_%s_group%d_THR%d_%s.npz' % (outdir, featureName, theTis, THR, FMfn))
    except:
	print('TFBS for %s  does not exist' % theTis)
	sys.exit()
    TF_allSNP = TF_allsnp_genes['TFBS_snps'][()]
    TF_allGenes = TF_allsnp_genes['TFBS_genes'][()]


    try:
	tp = pd.read_csv('%s/0815_noRestriction/TF_enrichment_%s_THR%d%s%s_TFBS_%s_torandom_idx%d.txt' % (outdir, featureName, THR, randomPattern, theTis, LMfn, idx), sep='\t')
    except:
    	outLier = readin_outliers(theTis)
        ts_pair = outLier[theTis]
        background_pair = outLier[-100]

	#[ts_pair, background_pair] = separate_pairs(ts_pair, bg_pair, countAs='both')
	result = collect_results(ts_pair, background_pair)
	result.to_csv('%s/0815_noRestriction/TF_enrichment_%s_THR%d%s%s_TFBS_%s_torandom_idx%d.txt' % (outdir, featureName, THR, randomPattern, theTis, LMfn, idx), sep='\t')




