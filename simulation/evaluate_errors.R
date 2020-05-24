evaluate_error <- function(trueL, trueF, lp, fp, sig_hits){
  if(dim(fp)[2] < dim(trueF)[2]){
    dif = (dim(trueF)[2] - dim(fp)[2])
    fp = cbind(fp, matrix(rep(0, dim(fp)[1] * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, dim(lp)[1] * dif), ncol = dif))
  }
  rankK = dim(trueF)[2]
  suppressWarnings(library('combinat'))
  ordering = permn(seq(1,rankK))
  f_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    f_cor[ord] = cor(as.vector(trueF), as.vector(fp[,ordering[[ord]]]))
  }
  
  l_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    l_cor[ord] = cor(as.vector(trueL), as.vector(lp[,ordering[[ord]]]))
  }
  
  ord_sum = f_cor + l_cor
  ord = which.max(abs(ord_sum))
  lp = lp[,ordering[[ord]]]
  fp = fp[,ordering[[ord]]]
  sig_hits = sig_hits[, ordering[[ord]]]
  
  lp = lp / matrix(rep(apply(lp, 2, function(x) max(abs(x))), dim(lp)[1]), nrow = dim(lp)[1], byrow = T)
  fp = fp / matrix(rep(apply(fp, 2, function(x) max(abs(x))), dim(fp)[1]), nrow = dim(fp)[1], byrow = T)
  colnames(fp) = seq(1,dim(fp)[2])
  rownames(fp) = seq(1, dim(fp)[1])
  
  fp[is.na(fp)] = 0 
  lp[is.na(lp)] = 0
  
  l_corr = cor(as.vector(trueL), as.vector(lp))
  f_corr = cor(as.vector(trueF), as.vector(fp))
  yhat = as.vector(abs(sig_hits) > 0)
  y = as.vector(abs(trueL) > 0)
  accuracy = sum(y == yhat) / length(y)

  precision = sum((y == yhat) * (yhat == 1)) / sum(yhat == 1)
  recall =  sum((y == yhat) * (yhat == 1)) / sum(y == 1)
  
  return(list(c(l_corr, f_corr, accuracy, precision, recall), lp, fp))
}

