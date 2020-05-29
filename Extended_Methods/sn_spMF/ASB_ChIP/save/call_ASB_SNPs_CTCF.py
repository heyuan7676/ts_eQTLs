import pandas as pd
import numpy as np
import os
import sys
import pdb
from scipy.stats import binom_test
from statsmodels.stats import multitest

from GLOBAL_VAR import *



alignmetn_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output'
SNP_in_TFBS_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output_GTExSNPs/'
outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/ChIP_ASB'



def get_ASB_ratio(fn, reads_filter = 10):
	print(fn)
	ChIP_seq_df = pd.read_csv('%s/%s' % (alignmetn_dir,fn), sep='\t', low_memory=False)
	reads_count = pd.DataFrame(ChIP_seq_df.groupby(['SNP_name', 'which_allele']).size())
	reads_count.columns = ['reads_count']

	reads_count['TFBS_SNP'] = [x[1] for x in reads_count.index]
	reads_count['SNP'] = [x[0] for x in reads_count.index]

	## exclude variants with reads < 10
	SNP_reads_count = reads_count.groupby('SNP').sum()
	print('Number of variants mapped to reads: ', len(SNP_reads_count))

	SNP_reads_count = SNP_reads_count[SNP_reads_count['reads_count'] > reads_filter] 
	print('Number of variants with > %d reads: ' % reads_filter, len(SNP_reads_count))

	reads_count = reads_count.loc[SNP_reads_count.index]

	## exclude reads on X chromosome
	reads_count = reads_count.iloc[np.where([not x.startswith('chrX') for x in reads_count['SNP']])[0]]

	### reads that map to variants tested in GTEx
	SNP_in_TFBS = pd.read_csv('%s/SNP_inTFBS_inGTEx_%s' % (SNP_in_TFBS_dir, fn), sep=' ', header=None)
	reads_count = reads_count.loc[SNP_in_TFBS[0]]
	print('Number of GTEx variants with mapped reads: ', len(set(reads_count['SNP'])))

	### exclude the locations that have only one read mapped to it
	#reads_count = reads_count[reads_count['reads_count'] > 1]

	### restrict to reads that map to both alleles
	variants_mapped_to = pd.DataFrame(reads_count.groupby('SNP').size())
	reads_count = reads_count.merge(variants_mapped_to, left_on='SNP', right_index=True)
	reads_count = reads_count[reads_count[0] == 2]
	print("Number of variants on heterozygous sites: ", len(reads_count) / 2)

        if restrict_to_Peaks:
                reads_count = restrict_to_peaks(fn, reads_count)

	total_mapped_reads = pd.DataFrame(reads_count.groupby('SNP')['reads_count'].sum())
	asb_permuted = []
	for i in range(len(total_mapped_reads)):
    		allele1 = np.sum([np.random.random() < 0.5 for k in range(total_mapped_reads.iloc[i]['reads_count'])])
    		asb_permuted.append(binom_test(float(allele1), total_mapped_reads.iloc[i]['reads_count']))

	all_asb = reads_count.groupby('SNP').apply(lambda x: binom_test(x['reads_count'].iloc[0], np.sum(x['reads_count'])))

	return [all_asb, asb_permuted, reads_count]



def restrict_to_peaks(fn, reads_count):
        ### restrict to reads in peaks
        mac_outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks/'
        peaks_dat = pd.read_csv('%s/%s_peaks.narrowPeak' % (mac_outdir, fn.replace('_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_formatted', '')), sep='\t', header =None, low_memory = False)
        for ri in range(1,23):
                peaks_dat.loc[peaks_dat[0] == str(ri),0] = ri
        reads_count['chr'] = [int(x.split('_')[0].replace('chr','')) for x in reads_count['SNP']]
        reads_count['location'] = [int(x.split('_')[1].replace('chr','')) for x in reads_count['SNP']]
        chromosome_length = reads_count.groupby('chr')['location'].max()
        chromosomes = {}
        for ri in range(1,23):
                chromosomes[ri] = np.zeros(chromosome_length[ri]+1)
                peaks_chr = peaks_dat[peaks_dat[0] == ri]
                for i in np.array(peaks_chr[[1,2]]):
                        chromosomes[ri][i[0]:i[1]+1] = 1

        inPeak = [chromosomes[int(x.split('_')[0].replace('chr',''))][int(x.split('_')[1])] for x in reads_count['SNP']]
        reads_count_inPeaks = reads_count.iloc[np.where(inPeak)[0]]
        print("Number of variants within peaks called by MACS: ", len(reads_count_inPeaks) / 2)
	return reads_count_inPeaks


from scipy.stats import chi2_contingency
def compare_to_random(group, asb_variants):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_ASB = len(np.intersect1d(tp[1], asb_variants))
    liver_total = float(len(np.intersect1d(tp[1], testable_SNPs))) - liver_ASB

    tp = pd.read_csv('%s/%s_outlierPairs_random_matched_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_random_SNPs = np.unique(tp[1])

    liver_random_ASB = float(len(np.intersect1d(liver_random_SNPs, asb_variants)))
    liver_random_total = float(len(np.intersect1d(liver_random_SNPs, testable_SNPs))) - liver_random_ASB

    if liver_total == 0:
        liver_total += 1
        liver_ASB += 1
    if liver_random_total == 0:
	liver_random_total += 1
        liver_random_ASB += 1

    print("Compare to random pairs: ", [[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])
    if liver_random_ASB == 0:
	OR = 1
	pv = 1
    else:
    	OR = (liver_ASB/liver_total / (liver_random_ASB/liver_random_total))
    	pv = chi2_contingency([[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])[1]
    return [liver_ASB, liver_total, liver_random_ASB, liver_random_total, OR, pv]




def compare_to_random_v2(group, variants_matched, asb_variants):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_ASB = len(np.intersect1d(tp[1], asb_variants))
    liver_total = float(len(np.intersect1d(tp[1], testable_SNPs))) - liver_ASB

    liver_random_SNPs = variants_matched

    liver_random_ASB = float(len(np.intersect1d(liver_random_SNPs, asb_variants)))
    liver_random_total = float(len(np.intersect1d(liver_random_SNPs, testable_SNPs))) - liver_random_ASB

    if liver_total == 0:
        liver_total += 1
        liver_ASB += 1
    if liver_random_total == 0:
        liver_random_total += 1
        liver_random_ASB += 1

    print("Compare to random pairs: ", [[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])
    if liver_random_ASB == 0:
        OR = 1
        pv = 1
    else:
        OR = (liver_ASB/liver_total / (liver_random_ASB/liver_random_total))
        pv = chi2_contingency([[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])[1]
    return [liver_ASB, liver_total, liver_random_ASB, liver_random_total, OR, pv]




def compared_to_background(group, asb_variants):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_ASB = len(np.intersect1d(tp[1], asb_variants))
    liver_total = float(len(np.intersect1d(tp[1], testable_SNPs))) - liver_ASB

    C = len(asb_variants)-liver_ASB
    D =  float(len(testable_SNPs) - C)

    print("Compare to background:" , [[liver_ASB, liver_total], [C, D]])
    if liver_total == 0:
	liver_total += 1
	liver_ASB += 1
    if D == 0:
	D += 1
	C += 1

    if C == 0:
	OR = 1
	pv = 1
    else:
    	OR = (liver_ASB/liver_total / (C/D))
    	pv = chi2_contingency([[liver_ASB, liver_total], [C, D]])[1]
    return [liver_ASB, liver_total, C,D, OR, pv]



def match_reads(reads_count_total):
    random_variant_groups = {}

    for group in range(23):
        tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
        factor_total = np.unique(np.intersect1d(tp[1], np.array(reads_count_total.index)))
        factor_reads_pool = reads_count_total.loc[factor_total]

    	variants_pool = np.array(list(set(np.array(reads_count_total.index)) - set(factor_total)))
    	variants_pool_asb = reads_count_total.loc[variants_pool]

    	random_variant = []

	random_fold = 10
    	for readsi in np.unique(factor_reads_pool):
		variant_to_match_number = np.sum(factor_reads_pool == readsi) * random_fold
        	pool_i = variants_pool_asb[np.abs(variants_pool_asb - readsi ) < 5]
        	if len(pool_i) < variant_to_match_number:
            		pool_i = variants_pool_asb[np.abs(variants_pool_asb - readsi) < 10]
        	if len(pool_i) < variant_to_match_number:
            		pool_i = variants_pool_asb[np.abs(variants_pool_asb - readsi) < 20]
        	if len(pool_i) < variant_to_match_number:
            		continue
        	pool_i = np.array(pool_i.index)
        	random_variant.append(np.random.choice(pool_i, variant_to_match_number, replace =False))

    	random_variant = np.unique([a for b in random_variant for a in b]  )
	random_variant_groups[group] = random_variant
    return random_variant_groups




if __name__ == '__main__':
	ChIP_type = sys.argv[1]
	reads_filter = int(sys.argv[2])
	restrict_to_Peaks = int(sys.argv[3])

	files = [x for x in os.listdir('/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output/') if ((ChIP_type in x) and x.endswith('formatted'))]
	variants, variants_matched  = [], []
	ASB, ASB_rd = [], []
	variants_matched_factors_files = {}	

	for fn in files:
		[abs_fn, abs_rd, reads_count]  = get_ASB_ratio(fn, reads_filter)
		reads_count_total = reads_count.groupby('SNP')['reads_count'].sum()
		## get random variants matched for number of reads
		tp = match_reads(reads_count_total)
		for key in tp.keys():
			if key not in variants_matched_factors_files.keys():
				variants_matched_factors_files[key] = list(tp[key])
			else:
				variants_matched_factors_files[key] = variants_matched_factors_files[key] + list(tp[key])
				variants_matched_factors_files[key] = np.unique(variants_matched_factors_files[key])

        	variants.append(np.array(abs_fn.index))
        	ASB.append(np.array(abs_fn))
		ASB_rd.append(abs_rd)


	variants = np.array([a for b in variants for a in b])
	ASB      = np.array([a for b in ASB      for a in b])
	ASB_rd   = np.array([a for b in ASB_rd   for a in b])
	testable_SNPs = np.unique(variants)

	df = pd.DataFrame({"SNP": variants, "ASB": ASB})
	df = pd.DataFrame(df.groupby('SNP')['ASB'].min())

        # empirical FDR
        p_thr = 0.05
	fdr = 1
        while fdr > 0.1:
                p_thr = p_thr * 0.9
                fdr = float(np.sum(ASB_rd < p_thr) ) / np.sum(ASB < p_thr)
	p_thr_01 = p_thr
        abs_fdr_01 = df.index[np.where(df['ASB'] < p_thr)]

        while fdr > 0.05:
                p_thr = p_thr * 0.9
                fdr = float(np.sum(ASB_rd < p_thr) ) / np.sum(ASB < p_thr)
	p_thr_005 = p_thr
        abs_fdr_005 = df.index[np.where(df['ASB'] < p_thr)]

	corrected_pv = multitest.multipletests(df['ASB'], method = 'fdr_bh')[1]
	abs_01  = df.index[np.where(corrected_pv < 0.1)[0]]
	abs_005 = df.index[np.where(corrected_pv < 0.05)[0]]


	print('Threshold for p-value to reach empirical FDR < 0.05 is %.4f' %  p_thr_01)
	print('Threshold for p-value to reach empirical FDR < 0.05 is %.4f' %  p_thr_005)

	print('Threshold for p-value to reach BH corrected FDR < 0.05 is %.4f' %  max(df.loc[abs_005]['ASB']))
	print('Threshold for p-value to reach BH corrected FDR < 0.1 is %.4f' %  max(df.loc[abs_01]['ASB']))
	
	print("Number of variants tesable in total:", len(df))
	print("Number of variants with ASB:", len(abs_01), len(abs_005), len(abs_fdr_01), len(abs_fdr_005))


	A,B,C,D = [], [], [], []
	OR, PV = [] , []
	group = []
	factor = []
	for g in list(range(23)):
		print(g)
		#print("Compare to random pairs:")
		for abs_variants_list in [abs_fdr_01, abs_fdr_005, abs_01, abs_005]:
			[a,b,c,d,orr,pvv] = compare_to_random(g, abs_variants_list)
			A.append(a)
			B.append(b)
			C.append(c)
			D.append(d)
    			OR.append(orr)
    			PV.append(pvv)


		## compare to the background
                #for abs_variants_list in [abs_fdr_01, abs_fdr_005, abs_01, abs_005]:
                #        [a,b,c,d,orr,pvv] = compared_to_background(g, abs_variants_list)
                #        A.append(a)
                #        B.append(b)
                #        C.append(c)
                #        D.append(d)
                #        OR.append(orr)
                #        PV.append(pvv)

                ## compare to the random pairs matched for read counts
                for abs_variants_list in [abs_fdr_01, abs_fdr_005, abs_01, abs_005]:
                        [a,b,c,d,orr,pvv] = compare_to_random_v2(g, variants_matched_factors_files[g], abs_variants_list)
                        A.append(a)
                        B.append(b)
                        C.append(c)
                        D.append(d)
                        OR.append(orr)
                        PV.append(pvv)


		group.append(["toRandom_fdr_01", "toRandom_fdr_005","toRandom_01", "toRandom_005", "toBG_fdr_01","toBG_fdr_005","toBG_01", "toBG_005"])
		factor.append([g] * 8)

	group = [a for b in group for a in b]
	factor = [a for b in factor for a in b]


	df = pd.DataFrame({"A":A, "B":B, "C":C, "D":D, "OR": OR, "pV": PV, "group": group, "factor": factor})


	if restrict_to_Peaks:
		outfn_pattern = '%s_Filtered%d_WithinPeaks' % (ChIP_type, reads_filter)
	else:
		outfn_pattern = '%s_Filtered%d_noMACSPeak'  % (ChIP_type, reads_filter)

	outfn_pattern = '%s_merged'% outfn_pattern

	outfn = '%s/Enrichment_analysis/ASB_%s_Enrichment.txt' % (outdir, outfn_pattern)
	print(outfn)
	df.to_csv(outfn, sep = '\t')


	fn = open('%s/ASB_SNPS/ASB_SNPs_%s_FDR01.txt' % (outdir, outfn_pattern ), 'w')
	for snpi in abs_01:
		fn.write('%s\n' % snpi)
	fn.close()

	fn = open('%s/ASB_SNPS/ASB_SNPs_%s_FDR005.txt' % (outdir, outfn_pattern), 'w')
        for snpi in abs_005:
                fn.write('%s\n' % snpi)
        fn.close()

        fn = open('%s/ASB_SNPS/ASB_SNPs_%s_emp_FDR01.txt' % (outdir, outfn_pattern ), 'w')
        for snpi in abs_fdr_01:
                fn.write('%s\n' % snpi)
        fn.close()

        fn = open('%s/ASB_SNPS/ASB_SNPs_%s_emp_FDR005.txt' % (outdir, outfn_pattern ), 'w')
        for snpi in abs_fdr_005:
                fn.write('%s\n' % snpi)
        fn.close()






