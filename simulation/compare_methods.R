suppressWarnings(library('plyr'))
suppressWarnings(library('dplyr'))
suppressWarnings(library('RColorBrewer'))
suppressWarnings(library('pheatmap'))
suppressWarnings(library('cowplot'))
suppressWarnings(library('reshape2'))
suppressWarnings(library('R.matlab'))
suppressWarnings(library('ggplot2'))
source("sn_spMF/readIn.R")
source('simulation/perform_MF_methods.R')
source('simulation/utils.R')


suppressWarnings(library(optparse))

option_list = list(make_option(c("-t", "--tau"), type = "double", default=1000, help="precision of the error"),
		   make_option(c("-s", "--seed"), type = "integer", default=1, help="seed of the random generator"))
opt = parse_args(OptionParser(option_list=option_list))



compare_methods <- function(tau, seed, K = 5, savedir = 'simulation/output', draw_result = F){
  dir.create(savedir, showWarnings = F)
  ## read in simulated input
  inputdir = 'simulation/input'
  file = paste0(inputdir, '/Input_tau', tau, '_seed', seed, '.RData')
  load(file)
  inputfnX = paste0(inputdir, '/Input_tau', tau, '_seed', seed, "_X.txt")
  inputfnW = paste0(inputdir, '/Input_tau', tau, '_seed', seed, "_W.txt")
  Data = readIn(inputfnX, inputfnW);
  X = Data[[1]]

  bnm = paste0('tau', tau, '_seed', seed)
  result = list()
  result[["NMF"]] = perform_MF_methods(X, l, f, method = 'NMF', rankK = K, bnm = bnm)
  result[["PCA"]] = perform_MF_methods(X, l, f, method = 'PCA', rankK = K, bnm = bnm)
  result[["softImpute"]] = perform_MF_methods(X, l, f, method = 'softImpute', rankK = K, bnm = bnm)
  result[["SSVD"]] = perform_MF_methods(X, l, f, method = 'SSVD', rankK = K, bnm = bnm)
  result[["PMA_cv1"]] = perform_MF_methods(X, l, f, method = 'PMA_cv1', rankK = K, bnm = bnm)
  result[["PMA_cv2"]] = perform_MF_methods(X, l, f, method = 'PMA_cv2', rankK = K, bnm = bnm)
  #result[["flashr_default"]] = perform_MF_methods(X, l, f, method = 'flashr_default', rankK = K, bnm = bnm)
  #result[["flashr_backfit"]] = perform_MF_methods(X, l, f, method = 'flashr_backfit', rankK = K, bnm = bnm)
  #result[["flashr_nn"]] = perform_MF_methods(X, l, f, method = 'flashr_nn', rankK = K, bnm = bnm)
  result[["sn_spMF"]] = perform_MF_methods(X, l, f, method = 'sn_spMF', rankK = K, bnm = bnm, readOnly = T)

  ## plot the learned factor matrices
  factors_matrices = NULL
  factors_matrices = rbind(format_f(f, 'True'))
  for(name in names(result)){
    factors_matrices = rbind(factors_matrices, format_f(result[[name]][[3]], name))
  }

  factors_matrices$method = factor(factors_matrices$method,
                                   levels=c("True", "sn_spMF","NMF",
                                            "flashr_nn",'flashr_backfit','flashr_default',
                                            'NBSFA', "softImpute", "PCA",
                                            "SSVD", "PMA_cv1", "PMA_cv2"))
  
  factors_matrices$dataPoint = factor(factors_matrices$dataPoint,
                                      levels = unique(factors_matrices$dataPoint))

  if(draw_result){  
  heatmap_g = ggplot(data = factors_matrices, aes(x = dataPoint, y = variable)) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = '#91CF60', mid = 'white', low = '#FC8D59', 
                         midpoint =  0, 
                         limits=c(min(factors_matrices$value),max(factors_matrices$value))) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    ylab("Factors") + xlab("Tissues")  + 
    theme(legend.position="top") + 
    labs(fill = "")  + 
    facet_wrap(~method, ncol = 3)
  
  #panel_width = unit(0.8,"npc") - sum(ggplotGrob(heatmap_g)[["widths"]][-3]) - unit(1,"line")
  #heatmap_g = heatmap_g + guides(fill= guide_colorbar(barwidth = panel_width))
  
  save_plot(paste0(savedir, "Learned_factor_matrices_tau", tau, '_seed',seed,'.png'), 
            heatmap_g, base_width = 6, base_height = 5)  
 
  }
 
  ### collect metrics
  factors_cors =  NULL
  for(name in names(result)){
    factors_cors = rbind(factors_cors, result[[name]][[1]])
  }
  
  factors_cors = as.data.frame(factors_cors)
  colnames(factors_cors) = c("l_corr", "f_corr", "accuracy", "precision", "recall")
  factors_cors$method = names(result)
  return(factors_cors)
}




compare_methods(opt$tau, opt$seed)
