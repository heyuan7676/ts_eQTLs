evaluate_error <- function(trueL, trueF, lp, fp, sig_hits){
  if(ncol(fp) < ncol(trueF)){
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
    sig_hits = cbind(sig_hits, matrix(rep(0, nrow(sig_hits) * dif), ncol = dif))
  }
  rankK = ncol(trueF)
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
  ord = which.max(ord_sum)
  lp = lp[,ordering[[ord]]]
  fp = fp[,ordering[[ord]]]
  sig_hits = sig_hits[, ordering[[ord]]]
  
  lp = lp / matrix(rep(apply(lp, 2, function(x) max(abs(x))), nrow(lp)), nrow = nrow(lp), byrow = T)
  fp = fp / matrix(rep(apply(fp, 2, function(x) max(abs(x))), nrow(fp)), nrow = nrow(fp), byrow = T)
  colnames(fp) = seq(1,ncol(fp))
  rownames(fp) = seq(1, nrow(fp))
  
  fp[is.na(fp)] = 0 
  lp[is.na(lp)] = 0
  
  l_corr = cor(as.vector(trueL), as.vector(lp))
  f_corr = cor(as.vector(trueF), as.vector(fp))
  
  ## u-eQTLs
  yhat = as.vector(abs(sig_hits[,1]) > 0)
  y = as.vector(abs(trueL[,1]) > 0)
  u_precision = sum((y == yhat) * (yhat == 1)) / sum(yhat == 1)
  u_recall =  sum((y == yhat) * (yhat == 1)) / sum(y == 1)

  ## u-eQTLs
  yhat = as.vector(abs(sig_hits[,seq(2, ncol(sig_hits))]) > 0)
  y = as.vector(abs(trueL[,seq(2, ncol(trueL))]) > 0)
  ts_precision = sum((y == yhat) * (yhat == 1)) / sum(yhat == 1)
  ts_recall =  sum((y == yhat) * (yhat == 1)) / sum(y == 1)
  
  return(list(c(l_corr, f_corr, u_precision, u_recall, ts_precision, ts_recall), lp, fp))
}

