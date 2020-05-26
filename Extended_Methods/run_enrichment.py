from collections import defaultdict,OrderedDict
import numpy as np
import sys
from scipy.stats import fisher_exact
import os
from gsea import readgenesets, gsea
import pandas as pd
import pdb

file_name = sys.argv[1] # File with test genes (one gene per line)
background_file_name = sys.argv[2] # File with background genes (one gene per line)
geneset_name = sys.argv[3] # filepath+filename of a geneset downloaded from msigdb. Need not parse msigdb file. Script takes care of parsing
save_name = sys.argv[4] # filepath+name for output file

""" Read test genes"""
fh = open(file_name,'r')
lines = fh.readlines()
fh.close()
test_genes = []
for line in lines:
	test_genes.append(line.strip('\n').split('\t')[1])


""" Read background genes"""
fh = open(background_file_name,'r')
lines = fh.readlines()
fh.close()
all_genes = []
for line in lines:
	all_genes.append(line.strip('\n').split('\t')[1])

		
all_genes = set(all_genes)
test_genes = set(test_genes)

[test_out, p, bf_p] = gsea(test_genes, all_genes, geneset_name)
order_idx = np.argsort(p)

bf_p = bf_p[order_idx]
test_out_ordered = OrderedDict()
for i in order_idx:
	key_i = test_out.keys()[i]
	test_out_ordered[key_i] = test_out[key_i]


""" Save result to save_name"""
f = open(save_name,'w')
#print >>f,"TestFile:", file_name, "\nBackgroundFile:", background_file_name, "\nGeneset:", geneset_name, "SaveFile:", save_name
print >>f,"\ttest_inset\ttest_notinset\tbackground_inset\tbackground_notinset\toddsratio\tpvalue\ttestgenes_in_set\tbonferroni-adjusted"
for k in range(len(bf_p)):
	i = test_out_ordered.keys()[k]
	print >>f, i,'\t','\t'.join(str(j) for j in test_out_ordered[i]), '\t', str(bf_p[k])
f.close()

#df = pd.read_csv(save_name, sep='\t', index_col = 0)
#df = df.sort_values('bonferroni-adjusted')
#print df.head()



