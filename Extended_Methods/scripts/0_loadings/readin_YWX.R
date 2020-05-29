library(data.table)
library(R.matlab)

readin_data <- function(FM_fn, prefix_input = "", Z_score = 0, startidx = "NONE"){

  ### read in data
  
  inputdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/'

  Y_filename = paste0(inputdir,  prefix_input , '_slope.txt')
  W_filename = paste0(inputdir,  prefix_input , '_se.txt')
  
  Y     = fread(Y_filename, sep='\t', header=T)
  Genes = Y$Gene
  SNPs  = Y$SNP
  keep_cols = colnames(Y)[3:51]
  Y         = Y[, ..keep_cols]
  Y         = data.frame(Y)
  
  W = fread(W_filename, sep='\t', header=T)
  W         = W[, ..keep_cols]
  W         = data.frame(W)

  
  W = 1/W
  
  if (startidx != "NONE"){
    startidx = as.numeric(startidx)
    gap = 100000
    #gap = 10
    idxS = startidx * gap + 1
    idxE = (startidx + 1) * gap
    if(idxS > dim(Y)[1]){
        quit() 
    }   
    if(idxE > dim(Y)[1]){
        idxE = dim(Y)[1]
    }
    Y = Y[idxS:idxE, ]
    W = W[idxS:idxE, ]

    Genes = Genes[idxS:idxE]
    SNPs  = SNPs[idxS:idxE]
  }
  
  
  N = dim(Y)[1]
  print(paste0(startidx, ': There are ', N, ' pairs tested'))
  
  
  ### read in the factor matrix
  if(FM_fn =='sn_spMF'){
  	fm_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph/'
  	fm_fn = paste0(fm_dir, 'SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110.mat')
  	fm = readMat(fm_fn)
  	X = fm$Xall[[1]][[1]]
  }else{
	X = read.table(paste0(FM_fn, '.txt'), sep='\t', header=T, stringsAsFactors = F)
  }

  X = as.matrix(X)
  X_max = apply(X, 2, function(x) max(abs(x)))
  X = t(t(X) / X_max)
  K = dim(X)[2]

  tissues = read.table('tissues.txt', sep='\t', header=F, stringsAsFactors = F)
  tissues = tissues$V1
  comp = apply(X, 2, function(x) tissues[which(abs(x)>0.01)])
  comp[[1]] = tissues

  if(Z_score == 1){
    Z = Y * W
    weighted_data = list()
    for(n in seq(1,N)){
      weighted_data[[n]] = list()
      weighted_data[[n]][[1]] = t(t(as.numeric(Z[n,])))
      weighted_data[[n]][[2]] = X
    }
  }else{
    ### construct data object
    weighted_data = list()
    for(n in seq(1,N)){
      weighted_data[[n]] = list()
      w = as.numeric(W[n, ])
      y = as.numeric(Y[n,])
      weighted_data[[n]][[1]] = t(t(y*w))
      weighted_data[[n]][[2]] = diag(w) %*% X
    }
  }
  return(list(Genes, SNPs, weighted_data, comp))
}




