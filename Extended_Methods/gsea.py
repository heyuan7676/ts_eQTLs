## Geneset Enrichment Analysis
## using gene symbol gmt files from  msigdb
## right tail Fisher's exact text

""" Read gene symbol gmt file from msigbdb"""
from statsmodels.stats import multitest
import numpy as np
import pdb

def readgenesets(filename):
	genesets={}
	with open(filename,'r') as fh:
		for line in fh:
			eachline = line.strip('\n').split('\t')
			#print(len(eachline[2::]))
			if(len(eachline[2::]) > 10 and len(eachline[2::]) < 1000):
				genesets[eachline[0]] = eachline[2::]
	return genesets

""" Fisher's exact test for enrichment"""
# Returns raw p-value
def gsea(testlist, background, gene_set):
	from scipy.stats import fisher_exact
	from collections import defaultdict,OrderedDict
	gene_sets = readgenesets(gene_set)
	size_gene_sets = len(gene_sets)
	enrich_out = OrderedDict()
	for i in gene_sets:
#		if len(gene_sets[i]) >= 50 and len(gene_sets[i]) <= 500:
		test_inset = len(set(testlist).intersection(gene_sets[i]))
		test_notinset = len(set(testlist)) - test_inset
		background_list = set(background) - set(testlist)
		#background_list = set(background)
		background_inset = len(set(background_list).intersection(gene_sets[i]))
		background_notinset = len(set(background_list)) - background_inset
		if test_inset + background_inset < 5:
			continue
		oddsratio, pvalue = fisher_exact([[test_inset,background_inset],[test_notinset,background_notinset]], alternative = 'greater')
		#oddsratio, pvalue = fisher_exact([[test_inset,background_inset],[test_notinset,background_notinset]])
		enrich_out[i] = tuple([test_inset, test_notinset, background_inset, background_notinset, oddsratio, pvalue, ','.join(set(testlist).intersection(gene_sets[i]))])
	print len(enrich_out)
	if len(enrich_out) == 0:
		return [[], []]
	pvalue = np.array([x[5] for x in enrich_out.values()])
	bf_adjusted = multitest.multipletests(pvalue,method = 'Bonferroni')[1]
	#bf_adjusted = np.array([min(1,p*70102) for p in pvalue])
	return [enrich_out, pvalue, bf_adjusted]
