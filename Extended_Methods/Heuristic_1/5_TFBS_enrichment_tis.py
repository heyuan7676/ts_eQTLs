import os
import sys
import numpy as np
import pandas as pd

from readin_X import readin_X
import pdb

from GLOBAL_VAR import *
from scipy.stats import fisher_exact
from scipy import stats
from scipy.stats import ranksums

def readin_outliers(theGROUP):
    ## outliers
    outLier = dict()
    outLier_df = pd.read_csv('%s/%s_outlierPairs_%sgroup%s.txt' % (pairdir, LMfn, randomPattern, theGROUP), sep='\t', header=None)
    outLier_df.columns = ['Gene', 'SNP']
    outLier[theGROUP] = outLier_df

    if background == 'shared':
    	outLier_df = pd.read_csv('%s/%s_outlierPairs_%sgroup%s.txt' % (pairdir, LMfn, randomPattern, 0), sep='\t', header=None)
    	outLier_df.columns = ['Gene', 'SNP']
    else:
        outLier_df = pd.read_csv('%s/%s_outlierPairs_%sgroup%s.txt' % (pairdir, LMfn, 'random_matched_', theGROUP), sep='\t', header=None)
        outLier_df.columns = ['Gene', 'SNP']


    # remove the ts-pairs
    #tsPairs = np.array(outLier[theGROUP].apply(lambda x: ':'.join((x['Gene'], x['SNP'])), axis=1))
    #backgroundPairs = np.array(outLier_df.apply(lambda x: ':'.join((x['Gene'], x['SNP'])), axis=1))
    #outLier_df.index = backgroundPairs
    #backgroundPaird_noTSeffect = set(backgroundPairs) - set(tsPairs)
    #outLier[rd_group] = outLier_df.loc[np.array(list(backgroundPaird_noTSeffect))]

    outLier[-1] = outLier_df

    return outLier




def match_SNPset_size(other_pairs, dist1, dist2):
    
        pair_pool = other_pairs.groupby('Gene')['SNP'].apply(list)
        matched_pairs = []
        for gi in dist2.index:
                loc_in_N = stats.percentileofscore(dist2, dist2[gi])
                N        = int(np.round(np.percentile(dist1, loc_in_N)))
                if N <= dist2[gi]:
                        matched_p = np.random.choice(pair_pool[gi], N, replace = False)
                else:
                        matched_p = np.random.choice(pair_pool[gi], N, replace = True)
                matched_pairs.append([':'.join((gi, s)) for s in matched_p])
                
        matched_pairs = [a for b in matched_pairs for a in b]
        other_pairs_new = pd.DataFrame([p.split(':') for p in matched_pairs])
        other_pairs = other_pairs_new.copy()
        other_pairs.columns = ['Gene', 'SNP']
        other_pairs = other_pairs.set_index('SNP', drop = False)
        dist2_new = other_pairs.groupby('Gene').size()
        
        return other_pairs




def separate_pairs(outLier_pairs, group, countAs = 'NONE'):

        ts_pair = outLier_pairs[group]
        background_pair = outLier_pairs[-1]     
        #print(len(set(ts_pair['Gene'])), len(ts_genes), len(set(background_pair['Gene'])), len(background_genes))
        
        # ts-pairs 
        if countAs == 'NONE':
	    ts_genes = np.array(list(set(ts_pair['Gene']) - set(background_pair['Gene'])))
            ts_pair = ts_pair.set_index('Gene', drop = False).loc[ts_genes].set_index('SNP', drop=False)

	    background_genes = np.array(list(set(background_pair['Gene']) - set(ts_pair['Gene'])))
            background_pair = background_pair.set_index('Gene', drop = False).loc[background_genes].set_index('SNP', drop=False)
            
        elif countAs == 'ts':
	    background_genes = np.array(list(set(background_pair['Gene']) - set(ts_pair['Gene'])))
            background_pair = background_pair.set_index('Gene', drop = False).loc[background_genes].set_index('SNP', drop=False)
            
	dist1 = ts_pair.groupby('Gene').size()
        dist2 = background_pair.groupby('Gene').size()

	if ranksums(dist1, dist2)[1] < 0.05:
		if   np.median(dist1) > np.median(dist2):
	    		ts_pair     = match_SNPset_size(ts_pair, dist2, dist1)
		elif np.median(dist1) < np.median(dist2):
            		background_pair = match_SNPset_size(background_pair, dist1, dist2)

	dist1_new = ts_pair.groupby('Gene').size()
        dist2_new = background_pair.groupby('Gene').size()

	print(group, np.median(dist1), np.median(dist2), np.median(dist1_new), np.median(dist2_new))
 
        return [ts_pair, background_pair]




def overlap_genes_between_genesets(ts_pair, background_pair, TF_genes, TF_snps):

        ts_pair_geneTFBS = ts_pair.set_index('Gene', drop=False).loc[np.unique(np.intersect1d(TF_genes, ts_pair['Gene']))]
        theSNP_TFBS = np.intersect1d(ts_pair_geneTFBS['SNP'], TF_snps)
        ts_pair_theSNP_TFBS = ts_pair_geneTFBS.set_index('SNP', drop = False).loc[theSNP_TFBS]

        background_pair_geneTFBS = background_pair.set_index('Gene', drop=False).loc[np.unique(np.intersect1d(TF_genes, background_pair['Gene']))]
        theSNP_TFBS = np.intersect1d(background_pair_geneTFBS['SNP'], TF_snps)
        background_pair_theSNP_TFBS = background_pair_geneTFBS.set_index('SNP', drop = False).loc[theSNP_TFBS]

	return [ts_pair_theSNP_TFBS, background_pair_theSNP_TFBS]

 
def find_asscio_pairs_for_snps_in_TFBS(the_pair, TF_genes, TF_snps):
    
        the_pair_geneTFBS = the_pair.set_index('Gene', drop=False).loc[np.unique(np.intersect1d(TF_genes, the_pair['Gene']))]
        #print('There are in total %d genes with any SNP in TFBS' % len(set(the_pair_geneTFBS['Gene'])))
    
        theSNP_TFBS = np.intersect1d(the_pair_geneTFBS['SNP'], TF_snps)
        the_pair_theSNP_TFBS = the_pair_geneTFBS.set_index('SNP', drop = False).loc[theSNP_TFBS]
        
        the_gene_theSNP_TFBS = set(the_pair_theSNP_TFBS['Gene'])
        the_gene_theSNP_nonTFBS = set(the_pair_geneTFBS['Gene']) - the_gene_theSNP_TFBS
        
        return [len(the_gene_theSNP_TFBS), len(the_gene_theSNP_nonTFBS)]


def find_asscio_pairs_for_snps_in_TFBS_nores(the_pair, TF_snps):

        theSNP_TFBS = np.intersect1d(the_pair['SNP'], TF_snps)
        the_pair_theSNP_TFBS = the_pair.set_index('SNP', drop = False).loc[theSNP_TFBS]

        the_gene_theSNP_TFBS = set(the_pair_theSNP_TFBS['Gene'])
        the_gene_theSNP_nonTFBS = set(the_pair['Gene']) - the_gene_theSNP_TFBS

        return [len(the_gene_theSNP_TFBS), len(the_gene_theSNP_nonTFBS)]

    

def contig_test(a, b, c, d):
        [a,b,c,d] = [float(a), float(b), float(c), float(d)]
        if 1:
            rs = fisher_exact([[a,b], [c,d]], alternative = 'greater')
        else:
            rs = chi2_contingency([[a,b],[c,d]])
        if c == 0:
                c += 1
                d += 1
                a += 1
                b += 1
        return [(a/b)/(c/d), rs[1]]




if __name__ == '__main__':
    theGROUP = int(sys.argv[1])
    featureName = sys.argv[2]
    background = sys.argv[3]
    THR = int(sys.argv[4])
    randomPattern = ""
    [X, Comp_tissues, tissues] = readin_X(FMfn)

    ## genes with SNPs in TFBS
    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/active_regions_features'
    try:
	TF_allsnp_genes = np.load('%s/TFBS_allsnpgenes/TFBS_allsnpgenes_%s_group%d_THR%d_%s_Esophagus_Mucosa.npz' % (outdir, featureName, theGROUP, THR, FMfn))
	#TF_allsnp_genes = np.load('%s/TFBS_allsnpgenes/TFBS_allsnpgenes_%s_group%d_THR%d_%s.npz' % (outdir, featureName, theGROUP, THR, FMfn))
    except:
	print('TFBS for Loading %d does not exist' % theGROUP)
	sys.exit()
    TF_allSNP = TF_allsnp_genes['TFBS_snps'][()]
    TF_allGenes = TF_allsnp_genes['TFBS_genes'][()]

    ## outliers
    outLier = readin_outliers(theGROUP)
    [ts_pair, background_pair] = separate_pairs(outLier, theGROUP, countAs='both')

    for tfi in np.sort(list(TF_allSNP.keys()))[:1]:
        if len(set(TF_allSNP[tfi])) == 0:
                continue
    	[df1, df2] = overlap_genes_between_genesets(ts_pair, background_pair, TF_allGenes[tfi],TF_allSNP[tfi])
	print(len(set(np.intersect1d(df1['Gene'], df2['Gene']))))

    print(ranksums(np.array(background_pair.groupby('Gene').size()), np.array(ts_pair.groupby('Gene').size())))


    ## eQTL in enhancer 
    outlier_fc = pd.read_csv('%s/%s_outlierSNPs_%s_group%d.txt' % (sigSNPfeaturedir, LMfn, featureName.split('_')[0], theGROUP), sep='\t', index_col = 0)
    outlier_fc.columns = [x.replace('chr1_','') for x in outlier_fc.columns]
    existing_tissue = np.intersect1d(Comp_tissues[theGROUP], outlier_fc.columns)


    ts_feature = outlier_fc.loc[ts_pair['SNP']][existing_tissue]
    eQTL_in_feature = np.zeros(len(ts_feature))
    for tis in existing_tissue:
    	eQTL_in_feature = np.array(['_'.join(featureName.split('_')[1:]) in xx for xx in np.array(ts_feature[tis])]) * 1
    ts_pair = ts_pair.iloc[np.where(eQTL_in_feature)[0]]


    if background == 'shared':
	outlier_fc = pd.read_csv('%s/%s_outlierSNPs_%s_group%d.txt' % (sigSNPfeaturedir, LMfn, featureName.split('_')[0], 0), sep='\t', index_col = 0)
    else:
    	outlier_fc = pd.read_csv('%s/%s_outlierSNPs_random_matched_%s_group%d.txt' % (sigSNPfeaturedir, LMfn, featureName.split('_')[0], theGROUP), sep='\t', index_col = 0)
    outlier_fc.columns = [x.replace('chr1_','') for x in outlier_fc.columns]

    background_feature = outlier_fc.loc[background_pair['SNP']][existing_tissue]
    eQTL_in_feature = np.zeros(len(background_feature))
    for tis in existing_tissue:
        eQTL_in_feature = np.array(['_'.join(featureName.split('_')[1:]) in xx for xx in np.array(background_feature[tis])]) * 1
    background_pair = background_pair.iloc[np.where(eQTL_in_feature)[0]]

    OR1, PV1, OR2, PV2, TF = [], [], [], [], []
    N1, N2, N3, N4, N5, N6 = [], [], [], [], [], []
    for tfi in np.sort(list(TF_allSNP.keys())):
        if len(set(TF_allSNP[tfi])) == 0:
                continue
        #print('Genomewide, there are %d Genes with any SNP on TFBS for %s' % (len(set(TF_allGenes[tfi])), tfi))
        [A, B] = find_asscio_pairs_for_snps_in_TFBS(ts_pair, TF_allGenes[tfi],TF_allSNP[tfi])
        [C, D] = find_asscio_pairs_for_snps_in_TFBS(background_pair, TF_allGenes[tfi],TF_allSNP[tfi])

        tfi = tfi.split('_')[1]
        if (B==0) or (D == 0):
            continue
        if contig_test(A, B, C, D)[1] < 1e-5:
            print(Comp_tissues[theGROUP], tfi, A, B, C,D, contig_test(A, B, C, D))
            
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

    if len(result) > 0:
    	result.to_csv('%s/TF_enrichment_%s_THR%d_%sgroup%d_TFBS_%s_to%s_Esophagus_Mucosa.txt' % (outdir, featureName, THR, randomPattern, theGROUP, LMfn, background), sep='\t')
	#result.to_csv('%s/TF_enrichment_%s_THR%d_%sgroup%d_TFBS_%s_to%s_restricted.txt' % (outdir, featureName, THR, randomPattern, theGROUP, LMfn, background), sep='\t')



