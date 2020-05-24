suppressWarnings(library(data.table))
source('sn_spMF/readIn.R')


readin_X <- function(FM_fn, factor_output_dir){
        ### read in the factor matrix
        ## if there are multiple runs, use the one whose factors being the most independent
        fx_dir = file.path(factor_output_dir, FM_fn)
        fx_files = list.files(fx_dir)
        fx_files = fx_files[grep('RData', fx_files)]


        min_obj = Inf
        for(fx in fx_files){
                load(paste0(fx_dir, '/', fx))
		obj = objective[length(objective)]
                print(c(fx, obj))
                if(obj < min_obj){
                        FactorM_optimal = FactorM
                        min_obj = obj
                        optimal_fx = fx
                }
        }
        print(paste0('Use the optimal factor matrix: ', optimal_fx))


        ## scale the factors to have the same magnitude (same max abs value)
        scale_max = median(apply(abs(FactorM_optimal), 2, max))
        FactorM_optimal = apply(abs(FactorM_optimal), 2, function(x) scale_max / max(x) * x) * sign(FactorM_optimal)
        K = ncol(FactorM_optimal)

        ## obtain the representative tissues/features for each factor
        tissues = rownames(FactorM_optimal)
        compositions = apply(FactorM_optimal, 2, function(x) tissues[which(x>0.01)])

        return(list(FactorM_optimal, compositions, tissues, optimal_fx))
}




readin_YWX <- function(FM_fn, xfn, wfn, factor_output_dir, F_C = NULL, Z_score = F, startidx = "NONE"){

	### read in data
	Data = readIn(Xfn = xfn, Wfn = wfn);
  	X = Data[['X']];
  	W = Data[['W']];
	Genes = Data[['Genes']];
	SNPs = Data[['SNPs']];
  
  	### split the dataset to run in parallel to speed up  
  	if (startidx != "NONE"){
    		startidx = as.numeric(startidx)
    		gap = 100000
    		idxS = startidx * gap + 1
    		idxE = (startidx + 1) * gap
    		if(idxS > nrow(X)){
        		quit() 
    		}   
    		if(idxE > nrow(X)){
        		idxE = nrow(X)
    		}
    		Y = Y[idxS:idxE, ]
    		W = W[idxS:idxE, ]

    		Genes = Genes[idxS:idxE]
    		SNPs  = SNPs[idxS:idxE]
  	}
  
  
	N = nrow(X)
  	print(paste0('There are ', N, ' pairs tested'))
  
  	### read in the factor matrix
	if(is.null(F_C)){
		F_C = readin_X(FM_fn, factor_output_dir)
		FactorM_optimal = F_C[[1]]
		compositions = F_C[[2]]
		tissues = F_C[[3]]
	}else{
		FactorM_optimal = F_C
		rownames(FactorM_optimal) = colnames(X)
        	## scale the factors to have the same magnitude (same max abs value)
        	scale_max = median(apply(abs(FactorM_optimal), 2, max))
        	FactorM_optimal = apply(abs(FactorM_optimal), 2, function(x) scale_max / max(x) * x) * sign(FactorM_optimal)
        	K = ncol(FactorM_optimal)

        	## obtain the representative tissues/features for each factor
        	tissues = rownames(FactorM_optimal)
        	compositions = apply(FactorM_optimal, 2, function(x) tissues[which(x>0.01)])
	}

  	if(Z_score){
		## in situations where no weights are used
    		z_mat = X * W
    	return(list(Genes, SNPs, z_mat, FactorM_optimal, compositions))
  	}else{
    		## construct data object
    		weighted_data = list()
    		for(n in seq(1,N)){
      			weighted_data[[n]] = list()
      			w = as.numeric(W[n, ])
      			y = as.numeric(X[n,])
      			weighted_data[[n]][[1]] = t(t(y*w))
      			weighted_data[[n]][[2]] = diag(w) %*% FactorM_optimal
    		}
  		return(list(Genes, SNPs, weighted_data, compositions, tissues))
  	}

}



readin_XWBF <- function(FM_fn, xfn, wfn, mapping_output_dir, factor_output_dir, F_C = NULL){

        ### read in data
        Data = readIn(Xfn = xfn, Wfn = wfn);
        X = Data[['X']];
        W = Data[['W']];
        Genes = Data[['Genes']];
        SNPs = Data[['SNPs']];

        N = nrow(X)
        print(paste0('There are ', N, ' pairs tested'))

        ### read in the factor matrix
        ### read in the factor matrix
        if(is.null(F_C)){
                F_C = readin_X(FM_fn, factor_output_dir)
                FactorM_optimal = F_C[[1]]
        }else{
                FactorM_optimal = F_C
                rownames(FactorM_optimal) = colnames(X)
                ## scale the factors to have the same magnitude (same max abs value)
                scale_max = median(apply(abs(FactorM_optimal), 2, max))
                FactorM_optimal = apply(abs(FactorM_optimal), 2, function(x) scale_max / max(x) * x) * sign(FactorM_optimal)
                K = ncol(FactorM_optimal)
        }

        ## read in the mapped betas
        B = read.csv(paste0(mapping_output_dir, 'mapping_', FM_fn, '_Loadings_beta_alpha0.05.txt'), sep='\t')

        return(list(X, W, B, FactorM_optimal))

}


