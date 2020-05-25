## flashr       * 0.6-7     2020-01-10 [1] Github (stephenslab/flashr@ae44e7d)
## ashr           2.2-39    2020-01-10 [1] Github (stephens999/ashr@e8a7abc)
## ebnm         * 0.1-24    2020-01-06 [1] Github (stephenslab/ebnm@308cb8a) 
## mixsqp         0.2-2     2019-10-16 [1] CRAN (R 3.5.2)      

my_init_fn <- function(Y, K = 1) {
  ret = flashr:::udv_si(Y, K)
  pos_sum = sum(ret$v[ret$v > 0])
  neg_sum = -sum(ret$v[ret$v < 0])
  if (neg_sum > pos_sum) {
    return(list(u = -ret$u, d = ret$d, v = -ret$v))
  } else
    return(ret)
}
flash_pipeline = function(data, rankK,m1='normal', m2='+uniform', ...) {
  ## current state-of-the art
  ## suggested by Jason Willwerscheid
  ## cf: discussion section of
  ## https://willwerscheid.github.io/MASHvFLASH/MASHvFLASHnn2.html
  ebnm_fn = "ebnm_ash"
  ebnm_param = list(l = list(mixcompdist = m1,
                             optmethod = "mixSQP"),
                    f = list(mixcompdist = m2,
                             optmethod = "mixSQP"))
  ##
  fl_g <- flashr:::flash_greedy_workhorse(data,
                                          var_type = "constant",
                                          ebnm_fn = ebnm_fn,
                                          Kmax = rankK,
                                          ebnm_param = ebnm_param,
                                          init_fn = "my_init_fn",
                                          stopping_rule = "factors",
                                          tol = 1e-3,
                                          verbose_output = "odF")
  
  
  fl_b <- flashr:::flash_backfit_workhorse(data,
                                           f_init = fl_g,
                                           var_type = "constant",
                                           ebnm_fn = ebnm_fn,
                                           ebnm_param = ebnm_param,
                                           stopping_rule = "factors",
                                           tol = 1e-3,
                                           verbose_output = "odF")
  return(fl_b)
}




tune_cor_PMA <- function(X, rankK, su, sv){
  N = dim(X)[1]
  Tn = dim(X)[2]
  random_idx = arrayInd(sample(N*Tn,N*Tn, replace =F),dim(X))
  gap = floor(N / 5)
  mres = 0
  for(i in seq(1,5)){
    rd_idx = random_idx[gap*(i-1)+1:gap*i, ]
    X_masked = X
    X_masked[rd_idx] = NA
    pmas = PMD(as.matrix(X_masked), K=rankK, sumabs = NULL, sumabsu = su, sumabsv = sv)
    X_hat = pmas$u %*% diag(pmas$d) %*% t(pmas$v)
    mres  = mres + sum((X_hat - X)[rd_idx]**2, na.rm = T) / (length(rd_idx) - sum(is.na((X_hat - X)[rd_idx])))
  }
  return(mres / 5)
}





format_f <- function(fm, method){
  ### format the factor matrix for plotting in ggplot2
  fm = data.frame(fm)
  P = dim(fm)[1]
  colnames(fm) = seq(1, dim(fm)[2])
  fm$dataPoint = seq(1,P)
  fm_melt = reshape2::melt(fm, id.vars = 'dataPoint')
  fm_melt$method = method
  return(fm_melt)
}
