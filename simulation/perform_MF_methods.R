source('simulation/evaluate_errors.R')
source('simulation/tune_parameters_function.R')


perform_MF_methods <- function(X, true_l, true_f, method, rankK = NULL, bnm = NULL, readOnly = F){
  
  available_methods = c('NBSFA', 'NMF', 'PCA', 
                        'softImpute', 'SSVD', 'PMA_cv1', 'PMA_cv2', 
                        'flashr_default','flashr_backfit', 'flashr_nn', 'sn_spMF')
  
  if(!method %in% available_methods){
    print(paste0('method should be one of ', available_methods))
    return()
  }
  
  print(paste0('Perform ', method))
  
  ### NBSFA
  if(method == 'NBSFA') {
	  suppressWarnings(library(R.matlab))
	  nbsfa = readMat(paste0('simulation/nsfa-master/code/results/X_',bnm, '.mat'))
	  l_hat = nbsfa$L
	  f_hat = t(nbsfa$F)
  }
  
  ### NMF
  if(method == 'NMF'){
    suppressWarnings(library(NMF))
    nmfs = nmf(abs(X), rank = rankK)
    f_hat = t(coef(nmfs))
    l_hat = basis(nmfs)
  }

  ### PCA
  if(method == 'PCA'){
    pcas = prcomp(X, scale. = F, center = F, rank. = rankK)
    f_hat = pcas$rotation 
    l_hat = pcas$x
  }


  ### sn_spMF
  if(method == 'sn_spMF'){
    sn_spMF_dir = paste0('simulation/output/', bnm, '/')
    res = tune_parameters(sn_spMF_dir, readOnly = readOnly)
    f_hat = res[[1]]
    l_hat = res[[2]]
  }

  ### softImpute
  if(method == 'softImpute'){
    suppressWarnings(library(softImpute))
    l1 = 10
    softimputes = softImpute(as.matrix(X), rank.max = rankK, lambda = l1)
    l_hat = softimputes$u
    f_hat = softimputes$v
    while(dim(f_hat)[2] < rankK){
      softimputes = softImpute(as.matrix(X), rank.max = rankK, lambda = l1)
      l_hat = softimputes$u
      f_hat = softimputes$v
      l1 = l1-0.5
    }
  }

  
  ### SSVD
  if(method == 'SSVD'){
    suppressWarnings(library(ssvd))
    ssvds = ssvd(as.matrix(X), method = 'method', r = rankK)
    l_hat = ssvds$u
    f_hat = ssvds$v
  }

  
  ### PMA
  if(method == 'PMA_cv1'){
    suppressWarnings(library(PMA))
    x = as.matrix(X)
    f_hat = NULL
    l_hat = NULL
    for(k in seq(1,rankK)){
      cv.out <- PMD.cv(x, type="standard", sumabss=seq(0.45,1, len=20))
      pmas = PMD(as.matrix(x), sumabs = cv.out$bestsumabs, sumabsu=NULL, sumabsv=NULL) ## sumabs: sparsity of v, sumabsu: sparsity of u
      x = x - pmas$u %*% pmas$d %*% t(pmas$v)
      f_hat = cbind(f_hat, pmas$v)
      l_hat = cbind(l_hat, pmas$u)
    }
  }
  
  ### PMA - cv2
  if(method == 'PMA_cv2'){
    suppressWarnings(library(PMA))
    result = NULL
    for(sv in seq(1,sqrt(5))){
      for(su in seq(1,10)){
        cv = tune_cor_PMA(as.matrix(X), rankK, su, sv)
        result = rbind(result, c(sv, su, cv))
      }
    }
  
    ops = result[which.min(result[,3]),]
    sv = ops[1]
    su = ops[2]
    pmas = PMD(as.matrix(X), K=rankK, sumabs = NULL, sumabsu = su, sumabsv = sv)
    f_hat = pmas$v
    l_hat = pmas$u
  }
  
  
  ### flashr - default
  if(method == 'flashr_default'){
    suppressWarnings(library(flashr))
    suppressWarnings(library(ashr))
    fmodel = flash(as.matrix(X), rankK)
    l_hat = fmodel$ldf$l
    f_hat = fmodel$ldf$f 
  }

  
  ### flashr - backfit
  if(method == 'flashr_backfit'){
    suppressWarnings(library(flashr))
    suppressWarnings(library(ashr))
    fmodel = flash_pipeline(as.matrix(X), rankK, m1 = 'normal', m2 = 'normal')
    l_hat = fmodel$ldf$l
    f_hat = fmodel$ldf$f 
    
  }

  ### flashr - non_negative factors
  if(method == 'flashr_nn'){
    suppressWarnings(library(flashr))
    suppressWarnings(library(ashr))
    fmodel = flash_pipeline(as.matrix(X), rankK)
    f_hat = fmodel$ldf$f  
    l_hat = fmodel$ldf$l    
  }

  ## get the significant data points
  xfn = paste0("simulation/input/Input_",bnm,"_X.txt")
  wfn = paste0("simulation/input/Input_",bnm,"_W.txt")
  if(method == 'NMF'){
	  sig_dp = l_hat
  }else{
	  sig_dp = fit_significant_hits(as.matrix(f_hat),paste0(bnm, '_', method), xfn, wfn)
  }

  metrics = evaluate_error(true_l, true_f, l_hat, f_hat, sig_dp)
  cat('\n')
  return(metrics)
}

