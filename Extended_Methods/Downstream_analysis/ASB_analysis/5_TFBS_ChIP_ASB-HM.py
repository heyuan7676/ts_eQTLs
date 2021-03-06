import pandas as pd
import numpy as np
import os
import sys


from GLOBAL_VAR import *



ChIP_type = "TFBS"
#hm = sys.argv[1]
fn = sys.argv[1]
theGROUP = 15
print(fn)

alignmetn_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/%s_ChIP_seq/STAR_output' % ChIP_type
SNP_in_TFBS_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/%s_ChIP_seq/STAR_output_GTExSNPs/' % ChIP_type
outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/ChIP_ASB'


## reads from ChIP-seq

#ChIP_seq_df = pd.read_csv('%s/liver_%s_seqAligned.sortedByCoord.out.bam.formatted' % (alignmetn_dir,hm), sep='\t', low_memory=False)
ChIP_seq_df = pd.read_csv('%s/%s' % (alignmetn_dir,fn), sep='\t', low_memory=False)
reads_count = pd.DataFrame(ChIP_seq_df.groupby(['SNP_name', 'which_allele']).size())
reads_count.columns = ['reads_count']

reads_count['TFBS_SNP'] = [x[1] for x in reads_count.index]
reads_count['SNP'] = [x[0] for x in reads_count.index]



## exclude variants with reads < 10
SNP_reads_count = reads_count.groupby('SNP').sum()
print('Number of variants mapped to reads: ', len(SNP_reads_count))

SNP_reads_count = SNP_reads_count[SNP_reads_count['reads_count'] > 10] 
print('Number of variants with > 10 reads: ', len(SNP_reads_count))

reads_count = reads_count.loc[SNP_reads_count.index]

## exclude reads on X chromosome
reads_count = reads_count.iloc[np.where([not x.startswith('chrX') for x in reads_count['SNP']])[0]]



### reads that map to variants tested in GTEx
SNP_in_TFBS = pd.read_csv('%s/SNP_inTFBS_inGTEx_%s' % (SNP_in_TFBS_dir, fn), sep=' ', header=None)
reads_count = reads_count.loc[SNP_in_TFBS[0]]
print('Number of GTEx variants with mapped reads: ', len(set(reads_count['SNP'])))


### restrict to reads that map to both alleles
variants_mapped_to = pd.DataFrame(reads_count.groupby('SNP').size())
reads_count = reads_count.merge(variants_mapped_to, left_on='SNP', right_index=True)
reads_count = reads_count[reads_count[0] == 2]
testable_SNPs = np.unique(reads_count['SNP'])
print("Number of variants on heterozygous sites: ", len(reads_count) / 2)




from scipy.stats import binom_test 
from statsmodels.stats import multitest
all_asb = reads_count.groupby('SNP').apply(lambda x: binom_test(x['reads_count'].iloc[0], np.sum(x['reads_count'])))
asb_variants = all_asb.index[np.where(multitest.multipletests(all_asb, method = 'fdr_bh')[0])]
print("Number of variants with ASB:", len(asb_variants))



from scipy.stats import chi2_contingency
def compare_to_random(group):
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_ASB = len(np.intersect1d(tp[1], asb_variants))
    liver_total = float(len(np.intersect1d(tp[1], testable_SNPs)))


    tp = pd.read_csv('%s/%s_outlierPairs_random_matched_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    liver_random_ASB = float(len(np.intersect1d(tp[1], asb_variants)))
    liver_random_total = float(len(np.intersect1d(tp[1], testable_SNPs)))

    print(group)
    print([[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])
    OR = (liver_ASB/liver_total / (liver_random_ASB/liver_random_total))
    pv = chi2_contingency([[liver_ASB, liver_total], [liver_random_ASB, liver_random_total]])[1]
    return [OR, pv]
    
    
OR, PV = [] , []
print("Compare to random pairs:")    
for g in [15]:
    [a,b] = compare_to_random(g)
    OR.append(a)
    PV.append(b)


df = pd.DataFrame({"OR": OR, "pV": PV})
outfn = '%s/ASB_%s_%s_toRandom.txt' % (outdir, ChIP_type, fn)
df.to_csv(outfn, sep = '\t')


OR = []
pv = []

tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, theGROUP), sep='\t', header = None)
liver_ASB = len(np.intersect1d(tp[1], asb_variants))
liver_total = float(len(np.intersect1d(tp[1], testable_SNPs)))


for group in range(0):
    if group == theGROUP:
        continue
        
    tp = pd.read_csv('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group), sep='\t', header = None)
    shared_ASB = len(np.intersect1d(tp[1], asb_variants))
    shared_total = float(len(np.intersect1d(tp[1], testable_SNPs)))

    #print([[liver_ASB, liver_total], [shared_ASB, shared_total]])
    OR.append((liver_ASB/liver_total / (shared_ASB/shared_total)))
    pv.append(chi2_contingency([[liver_ASB, liver_total], [shared_ASB, shared_total]])[1])
    
#idx = np.argsort(np.array(OR))
#OR = np.array(OR)
#df = pd.DataFrame({"OR": OR, "pv": pv})
#df.index = list(set(range(23)) - set([theGROUP]))

#outfn = '%s/ASB_%s.txt' % (outdir, fn)
#df.to_csv(outfn, sep = '\t')


#print('Enriched to eQTLs in %d other factors' % np.sum(multitest.multipletests(pv, method = 'fdr_bh')[0]))


