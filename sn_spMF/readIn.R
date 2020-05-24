#######################################################  
# Function to read in the data
# Please modify according to the data structure
####################################################### 

suppressWarnings(library(readr))

readIn <- function(Xfn, Wfn = NULL, rows_to_read = NULL){
	X = read_tsv(Xfn, col_type = cols());

	if(is.null(Wfn)){
		SE =  copy(X);
		SE[, seq(3,ncol(SE))] = 1;
	}else{
		SE = read_tsv(Wfn, col_type = cols());
	}

	if(sum(X$Gene != SE$Gene) > 0){
		stop("Error: Different genes in X and SE!")
	}

	if(sum(X$SNP != SE$SNP) > 0){
		stop("Error: Different SNPs in X and SE!")
	}

	Genes = X$Gene
	SNPs = X$SNP

	## remove the first two columns of gene names and SNP names
	X = X[, seq(3, ncol(X))];
	SE = SE[, seq(3, ncol(SE))];

	W = 1/SE;
	X[is.na(W)] = 0;
	W[is.na(W)] = 0;

	if(!is.null(rows_to_read)){
		idx = sample(seq(1, nrow(X)), rows_to_read, replace = F)
		X = X[idx, ]
		W = W[idx,]
	}

	return(list(X = X, W = W, Genes = Genes, SNPs = SNPs))
}





