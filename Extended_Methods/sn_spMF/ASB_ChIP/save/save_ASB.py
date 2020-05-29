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




if __name__ == '__main__':
	ChIP_type = sys.argv[1]
	reads_filter = int(sys.argv[2])
	restrict_to_Peaks = int(sys.argv[3])

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
	df.to_csv('%s/%s_merged_ASB.txt' % (outdir, outfn_pattern), sep='\t')

