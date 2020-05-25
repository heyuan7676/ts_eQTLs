source("sn_spMF/collect_results.R")
source("sn_spMF/readIn.R")
source("sn_spMF/fit_L.R")
source("simulation/fit_significant_hits.R")

tune_parameters <- function(outputdir, readOnly = F){
	opt = list()
	opt$outputdir = outputdir
	bnm = basename(opt$outputdir);

	if(readOnly){
		result = read.table(paste0('simulation/output/', bnm, '_choose_paras.txt'), sep='\t', header = T)
	}else{
		result = NULL
		for(K in c(5)){
			for(alpha1 in c(seq(1,10)/10, seq(1,10))){
				for(lambda1 in c(seq(1,10)/10, seq(1,10))){
					rowi = collect_results(opt$outputdir, K, alpha1, lambda1)
					if(!is.null(rowi)){
						result = rbind(result, rowi)
					}
				}
			}
		}

		result = as.data.frame(result)
		colnames(result) = c("K", "alpha1", "lambda1", "coph", "correlation", "nfactor", "run_optimal", "obj_optimal")
		result = result[order(result$coph), ]
		write.table(result, paste0('simulation/output/', bnm, '_choose_paras.txt'), sep='\t', quote = F, row.names = F)
	}


	res_filterd = find_optimal_Fm(result)
	res = res_filterd[[1]]
        K = res[1, "K"]
        a1 = res[1, "alpha1"]
        l1 = res[1, "lambda1"]
	run_idx = res[1, "run_optimal"]

	fn = res_filterd[[2]]
	ffn = paste0(opt$outputdir, '/sn_spMF_',fn,'/sn_spMF_FactorMatrix_',fn,'_Run',run_idx,'.txt')
	fM = read.table(ffn, sep='\t')


	# fit the loading matrix
        # read in the arguments
        option = list()
        option[['alpha1']]  = a1;
        option[['lambda1']] = l1;
        option[['K']]       = K;
        option[['iter']]  = 50;
        option[['convF']] = 0.02;
        option[['convO']] = 1;
        option[['disp']]  = F;

        bnm = basename(opt$outputdir);
        inputfnX = paste0("simulation/input/Input_",bnm,"_X.txt")
        inputfnW = paste0("simulation/input/Input_",bnm,"_W.txt")
        Data = readIn(inputfnX, inputfnW);

        lM = fit_L(Data[[1]], Data[[2]], as.matrix(fM), option);

	return(list(as.matrix(fM), as.matrix(lM)))
}



find_optimal_Fm <- function(result){
        result$nfactor = round(result$nfactor)
        res = result[result$nfactor == 5, ]
        if(sum(res$coph > 0.9) == 0){
                res = res[res$coph > 0.85,]
        }else{
                res = res[res$coph > 0.9, ]
        }
        res = res[order(res$correlation),]

        K = res[1, "K"]
        a1 = res[1, "alpha1"]

        l1 = res[1, "lambda1"]
		run_idx = res[1, "run_optimal"]

        fn = paste0('K', K, '_a1', a1, '_l1', l1)

	return(list(res, fn))
}



