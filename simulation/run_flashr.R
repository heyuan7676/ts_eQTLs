suppressWarnings(library('plyr'))
suppressWarnings(library('dplyr'))
suppressWarnings(library('reshape2'))
source("sn_spMF/readIn.R")
source('simulation/flashr_local.R')

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
  flashr_local(X, l, f, method = 'flashr_default', rankK = K, bnm = bnm)
  flashr_local(X, l, f, method = 'flashr_backfit', rankK = K, bnm = bnm)
  flashr_local(X, l, f, method = 'flashr_nn', rankK = K, bnm = bnm)
}


compare_methods(opt$tau, opt$seed) 
