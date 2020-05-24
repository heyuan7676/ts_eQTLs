catch_compensate_loadings_one_pair <- function(x, w, b, FactorM){

    idx = which(b!=0)
    if(length(idx) <= 1){ return(b)}

    idx = idx[idx!=1]
    compensate_loadings = c()
    N  = length(idx)

    for(i in seq(1, N)){
        k = idx[i]
	b_each = rep(0, length(b))
        b_each[k] = as.numeric(b[k])
        x_each = FactorM %*% as.matrix(b_each)
        x_each[is.na(x_each)] = 0

	compare_tissues = which(abs(x_each) > 0)
        x_each = x_each[compare_tissues]
        xw       = (x*w)[compare_tissues]

	compensate = rep(0, length(compare_tissues))
        for(j in seq(1, length(compensate))){
            if(sign(x_each[j]) * sign(xw[j]) == -1){
                compensate[j] = 1
		#print('Opposite signs')
	    }else if((sign(x_each[j]) * sign(xw[j]) == 1) & (abs(xw[j]) < 3)){
                compensate[j] = 1
		#print('Compensate effects')
	    }
	}

        if(sum(compensate==1) > length(compensate) * 0.8){
	    # for factor with multiple tissues, most of the tissues should have compensate effects to call the factor compensate
	    compensate_loadings = c(compensate_loadings, k)
	}

    if(length(compensate_loadings) > 0){
        b[compensate_loadings]  = 0
    }
    }

    return(b)

}



correct_B <- function(X, W, B, FactorM){

    B_correct = matrix(rep(0, nrow(B) * ncol(B)), ncol = ncol(B))
    for(ri in seq(1, nrow(X))){
        x = X[ri,]
        w = W[ri,]
        b = B[ri,]

        nonEmpty = which(b!=0)
        if(length(nonEmpty) == 0){next}

        ## re-fit using the non-zero coefficients
	nonnan = which(!is.na(x*w))
	Xw = diag(w[nonnan]) %*% FactorM[nonnan, ]
	lmmodel = lm(as.numeric((x*w)[nonnan]) ~ 0+Xw[, nonEmpty])
	fitted_b = rep(0, length(b))
	fitted_b[nonEmpty] = coef(lmmodel)

        b_correct = catch_compensate_loadings_one_pair(x, w, fitted_b, FactorM)
	B_correct[ri, ] = b_correct
    }

    B_correct = as.data.frame(B_correct)
    rownames(B_correct) = rownames(B)
    colnames(B_correct) = colnames(B)

    x1 = sum(B!=0) / nrow(B) / ncol(B)
    x2 = sum(B_correct!=0) / nrow(B) / ncol(B)
    print(paste0('Proportion of non-zero = ',x1, ' ; after correction = ', x2))

    return(B_correct)

}


