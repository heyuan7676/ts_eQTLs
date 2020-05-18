#######################################################  
# Function to read in the data
# Please modify according to the data structure
####################################################### 

readIn <- function(Xfn, Wfn = NULL, rows_to_read = NULL){
	X = read_tsv(Xfn);

	if(is.null(Wfn)){
		W = copy(X);
		W[, seq(3,ncol(W))] = 1;
	}else{
		W = read_tsv(Wfn);
	}

	if(sum(X$Gene != W$Gene) > 0){
		stop("Error: Different genes in X and W!")
	}

	if(sum(X$SNP != W$SNP) > 0){
		stop("Error: Different SNPs in X and W!")
	}

	## remove the first two columns of gene names and SNP names
	X = X[, seq(3, ncol(X))];
	W = W[, seq(3, ncol(W))];

	W = 1/W;
	X[is.na(W)] = 0;
	W[is.na(W)] = 0;

	if(!is.null(rows_to_read)){
		idx = sample(seq(1, nrow(X)), rows_to_read, replace = F)
		X = X[idx, ]
		W = W[idx,]
	}

	return(list(X = X, W = W))
}





