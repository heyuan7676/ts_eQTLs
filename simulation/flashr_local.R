source('simulation/evaluate_errors.R')
source('simulation/tune_parameters_function.R')
source('simulation/utils.R')

flashr_local <- function(X, true_l, true_f, method, rankK = NULL, bnm = NULL, readOnly = F){
  
  available_methods = c('flashr_default','flashr_backfit', 'flashr_nn')
  
  if(!method %in% available_methods){
    print(paste0('method should be one of ', available_methods))
    return()
  }
  
  print(paste0('Perform ', method))
  
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
  sig_dp  = fit_significant_hits(as.matrix(f_hat),paste0(bnm, '_', method), xfn, wfn)

  metrics = evaluate_error(true_l, true_f, l_hat, f_hat, sig_dp)
  save(metrics, file = paste0('flashr_local_results/', bnm, '_', method, '.RData'))
}

