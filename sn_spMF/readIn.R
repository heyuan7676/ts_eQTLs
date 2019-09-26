#######################################################  
# Function to read in the data
# X: Matrix of eQTL effect sizes 
# W: element wise inverse of the standard error
# Please modify according to the data structure
####################################################### 

readIn <- function(K, alpha1, lambda1, Xfn, Wfn){
	option = list()
	option[['alpha1']]  = alpha1;
	option[['lambda1']] = lambda1;
	option[['K']]       = K;

	option[['toobj']]   = 0.1;
	option[['disp']]    = TRUE;
	
	X = read.table(Xfn, sep='\t', header=T);
	W = read.table(Wfn, sep='\t', header=T);

	## remove the first two columns of gene names and SNP names
	#X = X[, seq(3, 51)];
	#W = W[, seq(3, 51)];

	W = 1/W;
	X[is.na(W)] = 0;
	W[is.na(W)] = 0;

	#set.seed(0)
	#idx = sample(seq(1, dim(X)[1]), 1000, replace = F)
	#X = X[idx, ]
	#W = 1/W[idx, ]
	#write.table(X, '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/test_data_X.txt', sep='\t')
	#write.table(W, '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/test_data_W.txt', sep='\t')

	Fn_basename = paste0('snspMF_K', K, '_a1', alpha1, '_l1', lambda1);
	return(list(X = X, W = W, option = option, Fn_basename = Fn_basename))
}





