datadir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/';
fn = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair';
Xfn = paste0(datadir, fn,'_slope.txt')
Wfn = paste0(datadir, fn,'_se.txt')
	
X = read.table(Xfn, sep='\t', header=T);
W = read.table(Wfn, sep='\t', header=T);

set.seed(0)
idx = sample(seq(1, nrow(X)), nrow(X), replace = F)
X = X[idx, ]
W = W[idx, ]

tissues = read.table('../data/tissues.txt', header = F, stringsAsFactors = F)
tissues = tissues$V1
colnames(X)[seq(3, ncol(X))] = tissues
colnames(W)[seq(3, ncol(W))] = tissues

write.table(X, '../data/test_data_X.txt', sep='\t', row.names = F, quote = F)
write.table(W, '../data/test_data_W.txt', sep='\t', row.names = F, quote = F)


