from GLOBAL_VAR import *
from METHOD_VAR import *
from collections import Counter

### read in gene annotation
gene_annotation = {}
gene_annotation_rev = {}
gene_annotation_fn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt','r')

for l in gene_annotation_fn.readlines():
        l = l.rstrip().split('\t')
        ensemID = l[0].replace('"','')
        geneName = l[5]
        gene_annotation[ensemID] = [geneName, l[1], l[2], l[3], l[6]]




shared_df = pd.read_csv('%s/%s_outlierPairs_group%s.txt' % (pairdir, LMfn, 0), sep='\t', header=None)
shared_df.columns = ['Gene', 'SNP']
shared_eGenes = np.unique(shared_df['Gene'])


all_eGenes = []
all_eGenes_dict = {}

for tis_index in range(1,23):
    outLier_df = pd.read_csv('%s/%s_outlierPairs_group%s.txt' % (pairdir, LMfn, tis_index), sep='\t', header=None)
    outLier_df.columns = ['Gene', 'SNP']
    all_eGenes.append(np.unique(outLier_df['Gene']))
    all_eGenes_dict[tis_index] = np.unique(outLier_df['Gene'])

all_eGenes = [a for b in all_eGenes for a in b]
gene_count = pd.DataFrame.from_dict(Counter(all_eGenes), orient='index')



for topGene_N in list(range(4,8)) + [30]:
    gene_set = gene_count[gene_count[0] < topGene_N]
    shared_gene_set = np.intersect1d(np.unique(gene_set.index), np.unique(shared_eGenes))
    fn = open('%s/%s_topGene%d_group%d.txt' % (pairdir, LMfn, topGene_N, 0), 'w')
    for g in shared_gene_set:
        fn.write('%s\t%s\n' % (g, gene_annotation[g][0]))
    fn.close()

    gene_set = list(set(gene_set.index) - set(shared_eGenes))
    
    x = []
    for tis_index in all_eGenes_dict.keys():
        fn = open('%s/%s_topGene%d_group%d.txt' % (pairdir, LMfn, topGene_N, tis_index), 'w')
        tpset = np.intersect1d(all_eGenes_dict[tis_index], gene_set)
        x.append(len(tpset))
        for g in tpset:
            fn.write('%s\t%s\n' % (g, gene_annotation[g][0]))
        fn.close()
        
    print(topGene_N, np.mean(x))



##### Tested genes


tested_gene_dir = '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL'

for tis_index in range(len(Comp_tissues)):
    print(tis_index)
    tp = pd.DataFrame()
    for tis in Comp_tissues[tis_index]:
        print(tis)
        tested_genes_dat = pd.read_csv('%s/%s.v8.egenes.txt' % (tested_gene_dir, tis), sep='\t')
    
        protein_coding_idx = np.where([gene_annotation[g][4] == 'protein_coding' for g in tested_genes_dat['gene_id']])[0]
        tested_genes_dat = tested_genes_dat.iloc[protein_coding_idx]
        tp.append(tested_genes_dat)

    tested_genes_dat[['gene_id', 'gene_name']].to_csv('%s/Revision/tested_genes/%s_genes_%s.txt' % (inputdatadir, tis_index, FMfn), sep='\t', index = False, header=None)
    
    
    
