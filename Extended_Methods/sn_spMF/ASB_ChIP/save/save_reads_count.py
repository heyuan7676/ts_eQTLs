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


def save_ASB_ratio(fn, reads_filter = 10):
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

        save_read_counts_fn = '%s/reads_count/%s_reads_filter_%d_peaks%d.txt' % (outdir, fn, reads_filter, restrict_to_Peaks)
        reads_count.to_csv(save_read_counts_fn, sep='\t')



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


if __name__ == '__main__':
	fn = sys.argv[1]
	reads_filter = int(sys.argv[2])
	restrict_to_Peaks = int(sys.argv[3])

	save_ASB_ratio(fn, reads_filter)



