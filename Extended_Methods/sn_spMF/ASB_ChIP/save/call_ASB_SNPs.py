import pandas as pd
import numpy as np
import os
import sys
import pdb
from scipy.stats import binom_test
from statsmodels.stats import multitest
from collections import Counter

from GLOBAL_VAR import *



alignmetn_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output'
SNP_in_TFBS_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output_GTExSNPs/'
outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/ChIP_ASB'



def get_ASB_ratio(save_read_counts_fn, reads_filter = 10):
	reads_count = pd.read_csv('%s/reads_count/%s' % (outdir, save_read_counts_fn), sep='\t', index_col = [0,1])
	
	total_mapped_reads = pd.DataFrame(reads_count.groupby('SNP')['reads_count'].sum())
	asb_permuted = []
	for i in range(len(total_mapped_reads)):
    		allele1 = np.sum([np.random.random() < 0.5 for k in range(total_mapped_reads.iloc[i]['reads_count'])])
    		asb_permuted.append(binom_test(float(allele1), total_mapped_reads.iloc[i]['reads_count']))

	all_asb = reads_count.groupby('SNP').apply(lambda x: binom_test(x['reads_count'].iloc[0], np.sum(x['reads_count'])))
	all_asb = pd.DataFrame(all_asb)

	reads_count_total = reads_count.groupby('SNP')['reads_count'].sum()	
	reads_count_total = pd.DataFrame(reads_count_total)

	return [all_asb, reads_count_total]



from scipy.stats import chi2_contingency

def compare_test(tis_ASB, tis_total, random_pairs, asb_variants):

    tis_random_ASB = float(len(np.intersect1d(random_pairs, asb_variants)))
    tis_random_total = float(len(np.intersect1d(random_pairs, testable_SNPs))) - tis_random_ASB

    if ((tis_total == 0) or (tis_random_total == 0)):
        tis_total += 1
        tis_ASB += 1
        tis_random_total += 1
        tis_random_ASB += 1

    print("Compare to random pairs: ", [[tis_ASB, tis_total], [tis_random_ASB, tis_random_total]])
    OR = (tis_ASB/tis_total / (tis_random_ASB/tis_random_total))
    pv = chi2_contingency([[tis_ASB, tis_total], [tis_random_ASB, tis_random_total]])[1]

    return [tis_ASB, tis_total, tis_random_ASB, tis_random_total, OR, pv]



def compare_to_random(group, asb_variants):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    ture_pairs = np.unique(tp[1])
    tis_ASB = len(np.intersect1d(ture_pairs, asb_variants))
    tis_total = float(len(np.intersect1d(ture_pairs, testable_SNPs))) - tis_ASB

    random_tp = pd.read_csv('%s/%s_outlierPairs_random_matched_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    random_pairs = np.unique(random_tp[1])
    R1 = compare_test(tis_ASB, tis_total, random_pairs, asb_variants)

    R2 = compare_test(tis_ASB, tis_total, matched_random_variants[group], asb_variants)

    return [R1, R2]



def compared_to_background(group, asb_variants):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    tis_ASB = len(np.intersect1d(tp[1], asb_variants))
    tis_total = float(len(np.intersect1d(tp[1], testable_SNPs))) - tis_ASB

    C = len(asb_variants)-tis_ASB
    D =  float(len(testable_SNPs) - C)

    print("Compare to background:" , [[tis_ASB, tis_total], [C, D]])
    if tis_total == 0:
	tis_total += 1
	tis_ASB += 1
    if D == 0:
	D += 1
	C += 1

    if C == 0:
	OR = 1
	pv = 1
    else:
    	OR = (tis_ASB/tis_total / (C/D))
    	pv = chi2_contingency([[tis_ASB, tis_total], [C, D]])[1]
    return [tis_ASB, tis_total, C,D, OR, pv]



def match_reads(reads_count_total):
    random_variant_groups = {}

    for group in groupList:
        tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
        factor_total = np.unique(np.intersect1d(np.unique(tp[1]), np.array(reads_count_total['SNP'])))
        factor_reads_pool = reads_count_total.loc[factor_total]

	print(group, Counter(factor_reads_pool['which_fn']))

	random_variant = []
	for reads_for_one_fn in factor_reads_pool.groupby('fn'):

		### info for the true variants
		one_fn = reads_for_one_fn[0]
		reads_for_one_fn = reads_for_one_fn[1]	

		### info for all testable variants -- used as pool for random variants
		reads_pool_for_this_fn = reads_count_total[reads_count_total['fn'] == one_fn]
    		variants_pool = np.array(list(set(np.array(reads_pool_for_this_fn.index)) - set(reads_for_one_fn.index)))
    		variants_pool_asb = reads_count_total.loc[variants_pool]

    		random_variant_fni = []

		random_fold = 2
    		for readsi in np.unique(reads_for_one_fn['total_count']):
			variant_to_match_number = np.sum(reads_for_one_fn['total_count'] == readsi) * random_fold

			logreadsi = np.log10(readsi)
        		pool_i = variants_pool_asb[np.asb(variants_pool_asb['log_total_count'] - logreadsi ) < 0.01]

        		if len(pool_i) < variant_to_match_number:
            			pool_i = variants_pool_asb[np.asb(variants_pool_asb['log_total_count'] - logreadsi) < 0.1]
        		if len(pool_i) < variant_to_match_number:
            			pool_i = variants_pool_asb[np.asb(variants_pool_asb['log_total_count'] - logreadsi) < 0.2]
			if len(pool_i) < variant_to_match_number:
				pool_i = variants_pool_asb[np.asb(variants_pool_asb['log_total_count'] - logreadsi) < 0.5]

			pool_i = np.array(pool_i.index)
        		if len(pool_i) < variant_to_match_number:
				rd_variants = np.random.choice(pool_i, variant_to_match_number, replace = True)
			else:
				rd_variants = np.random.choice(pool_i, variant_to_match_number, replace =False)	
        		random_variant_fni.append(rd_variants)

    		random_variant_fni = np.unique([a for b in random_variant_fni for a in b]  )
		random_variant.append(random_variant_fni)

	random_variant = [a for b in random_variant for a in b]	
	random_variant_groups[group] = random_variant

	#print(group, len(factor_total), len(random_variant))

    return random_variant_groups




if __name__ == '__main__':
	ChIP_type = sys.argv[1]
	reads_filter = int(sys.argv[2])
	restrict_to_Peaks = int(sys.argv[3])
	groupList = range(23)

	files = [x for x in os.listdir('%s/reads_count' % outdir) if ((ChIP_type in x) and ('reads_filter_%d' % reads_filter in x) and ('peaks%d' % restrict_to_Peaks in x))]
	variants, variants_matched  = [], []
	ASB, ASB_rd = [], []

	asb_fn = pd.DataFrame()
	reads_count = pd.DataFrame()

	for idx in range(len(files)):
		fn = files[idx]
		[asb_fn_i, reads_count_i]  = get_ASB_ratio(fn)
		reads_count_i['fn'] = idx
		asb_fn_i['fn'] = idx

		asb_fn = asb_fn.append(asb_fn_i)
		reads_count = reads_count.append(reads_count_i)

	asb_fn['SNP'] = asb_fn.index
	asb_fn.columns = ['ASB', 'fn', 'SNP']
	asb_fn['total_count'] = reads_count['reads_count']


	which_fn = asb_fn.groupby('SNP').apply(lambda x: x['fn'][np.argmin(np.array(x['ASB']))])
	which_fn = pd.DataFrame(which_fn)

	df = asb_fn.merge(which_fn, left_on = 'SNP', right_on = 'SNP')
	df.columns = ['ASB', 'fn', 'SNP', 'total_count', 'which_fn']
	df['log_total_count'] = np.log10(df['total_count'])
	df.index = df['SNP']
	df = df[df['fn'] == df['which_fn']]

	outfn_pattern = '%s_Filtered%d_peaks%d' % (ChIP_type, reads_filter, restrict_to_Peaks)
	df = pd.read_csv('%s/%s_merged_ASB.txt' % (outdir, outfn_pattern), sep='\t')


	testable_SNPs = np.unique(df.index)

	## get random variants matched for number of reads
	matched_random_variants = match_reads(df)

        # empirical FDR
        #p_thr = 0.05
	#fdr = 1
        #while fdr > 0.1:
        #        p_thr = p_thr * 0.9
        #        fdr = float(np.sum(ASB_rd < p_thr) ) / np.sum(ASB < p_thr)
	#p_thr_01 = p_thr
        #asb_fdr_01 = df.index[np.where(df['ASB'] < p_thr)]

        #while fdr > 0.05:
        #        p_thr = p_thr * 0.9
        #        fdr = float(np.sum(ASB_rd < p_thr) ) / np.sum(ASB < p_thr)
	#p_thr_005 = p_thr
        #asb_fdr_005 = df.index[np.where(df['ASB'] < p_thr)]

	corrected_pv = multitest.multipletests(df['ASB'], method = 'fdr_bh')[1]
	asb_01  = df.index[np.where(corrected_pv < 0.1)[0]]
	asb_005 = df.index[np.where(corrected_pv < 0.05)[0]]


	#print('Threshold for p-value to reach empirical FDR < 0.05 is %.4f' %  p_thr_01)
	#print('Threshold for p-value to reach empirical FDR < 0.05 is %.4f' %  p_thr_005)

	#print('Threshold for p-value to reach BH corrected FDR < 0.05 is %.4f' %  max(df.loc[asb_005]['ASB']))
	#print('Threshold for p-value to reach BH corrected FDR < 0.1 is %.4f' %  max(df.loc[asb_01]['ASB']))
	
	print("Number of variants tesable in total:", len(df))
	print("Number of variants with ASB:", len(asb_01), len(asb_005))


	group, factor, result_all = [], [], []
	for g in groupList:
		print(g)
		for asb_variants_list in [asb_01, asb_005]:
			[pvv1, pvv2] = compare_to_random(g, asb_variants_list)
			result_all.append(pvv1)
			result_all.append(pvv2)

		group.append(["matchedMAFTSS_01", "matchedRC_01","matchedMAFTSS_005", "matchedRC_005"])
		factor.append([g] * 4)

	group = [a for b in group for a in b]
	factor = [a for b in factor for a in b]
	result_all = pd.DataFrame(result_all)
	result_all.columns = ['A', 'B', 'C', 'D', 'OR', 'pV']
	result_all['group'] = group
	result_all['factor'] = factor

	outfn = '%s/Enrichment_analysis/ASB_%s_Enrichment.txt' % (outdir, outfn_pattern)
	print(outfn)
	result_all.to_csv(outfn, sep = '\t')

	fn = open('%s/ASB_SNPS/ASB_SNPs_%s_FDR01.txt' % (outdir, outfn_pattern ), 'w')
	for snpi in asb_01:
		fn.write('%s\n' % snpi)
	fn.close()

	fn = open('%s/ASB_SNPS/ASB_SNPs_%s_FDR005.txt' % (outdir, outfn_pattern), 'w')
        for snpi in asb_005:
                fn.write('%s\n' % snpi)
        fn.close()






