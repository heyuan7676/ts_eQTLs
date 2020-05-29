import pandas as pd
import numpy as np
import os
import sys


from GLOBAL_VAR import *


[X, Comp_tissues, tissues] = readin_X(FMfn)
alignmetn_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output'



def format_filtered_bam(fn):
    tp_sam_df = pd.read_csv('%s/%s' % ( alignmetn_dir, fn),  sep='\t', header=None, usecols=[2,17,18], low_memory=False)
    tp_sam_df.columns = ['chromosome', 'which_allele', 'location']

    location = [[int(xi) for xi in x.split(',')[1:]] for x in tp_sam_df['location']]
    which_allele = [[int(xi) for xi in x.split(',')[1:]] for x in tp_sam_df['which_allele']]

    chromosome = []
    for token in zip(np.array(tp_sam_df['chromosome']), [len(x) for x in location]):
        chromosome.append([token[0]] * token[1])
    
    chromosome = [a for b in chromosome for a in b]
    location   = [a+1 for b in location   for a in b]
    which_allele = [a for b in which_allele for a in b]


    snpname = []
    for token in zip(chromosome, location):
        snpname.append('chr%s_%d' % (token[0], token[1]))
    
    
    df = pd.DataFrame({"SNP_name": snpname, "chromosome": chromosome, 
                       "location": location, 'which_allele': which_allele})


    df.to_csv('%s/%s_formatted' % (alignmetn_dir, fn), sep='\t', index = False)



fn=sys.argv[1]
format_filtered_bam(fn)

