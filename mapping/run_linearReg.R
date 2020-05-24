compute_p <- function(dataPoint, tissues, compositions){
     #### run linear regression
     factor_exist = seq(1, dim(dataPoint[[2]])[2])
 
     # remove factors where no tissue has non-NA (0) values in it, and assign 0 to tissues with NA value
     if(sum((dataPoint[[1]] == 0)) > 0){
         valid_idx = which((dataPoint[[1]] != 0))
	 invalid_idx = which((dataPoint[[1]] == 0))
         data_exist = tissues[valid_idx]
         factor_exist = which(sapply(compositions, function(x) length(intersect(data_exist, x))>length(x)/2))

	 dataPoint[[1]][invalid_idx] = 0
	 if(dim(dataPoint[[1]])[1] != length(tissues)){
	 	dataPoint[[1]] = t(dataPoint[[1]])
	 }
         dataPoint[[2]] = dataPoint[[2]][, factor_exist]
	 if(length(dataPoint[[2]]) == length(tissues)){
		dataPoint[[2]] = t(t(dataPoint[[2]]))
	 }
	 dataPoint[[2]][invalid_idx,] = 0

	 # remove the shared factor if there are other factors become co-linear with it (because of NA values)
	 if(1 %in% factor_exist){
	 	temp_d2 = as.matrix(dataPoint[[2]][, 2:length(factor_exist)])
	 	if(sum(dataPoint[[2]][, 1] != 0) == max(apply(temp_d2, 2, function(x) sum(x!=0)))){
		 	dataPoint[[2]] = dataPoint[[2]][, seq(2, length(factor_exist))]
                 	factor_exist = factor_exist[factor_exist!=1]
         	}
	 }

         if(length(factor_exist) == 0){
                 beta = rep(0, length(compositions))
                 pv = rep(-1, length(compositions))
                 return(list(beta, pv))
         }
     }


     # remove the shared factor if the sign across available tissues disagree
     if(1 %in% factor_exist){
     	 if((sum(dataPoint[[1]] > 0, na.rm = T) > length(dataPoint[[1]])/3) & (sum(dataPoint[[1]] > 0, na.rm = T) < length(dataPoint[[1]]) / 3  * 2)){
		dataPoint[[2]] = dataPoint[[2]][, seq(2, length(factor_exist))]
         	factor_exist = factor_exist[factor_exist!=1]
	 }
     }


     lmfit = lm(dataPoint[[1]]~0+dataPoint[[2]], na.action = na.omit)
     beta = unlist(coef(lmfit))
     pv = coef(summary(lmfit))[,4]

     if(length(pv) < length(beta)){
	pv[names(beta)[which(is.na(beta))]] = -1
	names(pv) = matrix(unlist(strsplit(names(pv), "]]")), nrow=length(pv), byrow = T)[,2]
	pv = pv[order(as.numeric(names(pv)))]
     }

     beta[is.na(beta)] = 0
     pv[is.na(pv)] = -1
     if(length(pv) < length(compositions)){
         beta_formatted = rep(0, length(compositions))
         pv_formatted = rep(-1, length(compositions))
         beta_formatted[factor_exist] = beta
         pv_formatted[factor_exist] = pv
         beta = beta_formatted
         pv = pv_formatted
     }
     return(list(as.numeric(beta), as.numeric(pv)))
}




run_regression <- function(weighted_data, tissues, compositions){
    options(warn=-1)
    result = lapply(weighted_data, function(x) compute_p(x, tissues, compositions))
    pValues = c()
    Betas   = c()

    for(r in result){
        Betas = rbind(Betas, r[[1]])
        pValues = rbind(pValues, r[[2]])
    }

    Betas = as.data.frame(Betas)
    pValues = as.data.frame(pValues)

    return(list(Betas, pValues))

}


