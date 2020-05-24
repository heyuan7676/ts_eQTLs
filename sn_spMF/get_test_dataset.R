datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/';
fn = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair';
Xfn = paste0(datadir, fn,'_slope.txt')
Wfn = paste0(datadir, fn,'_se.txt')
	
X = read.table(Xfn, sep='\t', header=T);
W = read.table(Wfn, sep='\t', header=T);

tissues = read.table('data/tissues.txt', header = F, stringsAsFactors = F)
tissues = tissues$V1
colnames(X)[seq(3, ncol(X))] = tissues
colnames(W)[seq(3, ncol(W))] = tissues

set.seed(1)
tissues_for_demo = sort(sample(tissues, 20, replace = F))
X = X[, c("Gene", "SNP", tissues_for_demo)]
W = W[, c("Gene", "SNP", tissues_for_demo)]

write.table(X, 'data/test_data_X_all.txt', sep='\t', row.names = F, quote = F)
write.table(W, 'data/test_data_SE_all.txt', sep='\t', row.names = F, quote = F)


set.seed(0)
idx = sample(seq(1, nrow(X)), round(nrow(X)*0.1), replace = F)
X = X[idx, ]
W = W[idx, ]

write.table(X, 'data/test_data_X.txt', sep='\t', row.names = F, quote = F)
write.table(W, 'data/test_data_SE.txt', sep='\t', row.names = F, quote = F)


