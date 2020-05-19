
source("cophenet.R")
suppressWarnings(library(plyr))

collect_results <- function(K, alpha1, lambda1){

	# read in the arguments
	outputdir = '../output'
	Fn = paste0('K', K, '_a1', alpha1, '_l1', lambda1); 
	outputdir = file.path(outputdir, paste0('sn_spMF_learn_hp_snspMF_', Fn));

	F_all = list()
	fcorr_all = list()
	iii = 1
	coph_test = c()
	for(run_idx in seq(1,30)){
		outputFn = paste0(outputdir, '/sn_spMF_', Fn, '_Run',run_idx, '.RData');
		factorFn = paste0(outputdir, '/sn_spMF_FactorMatrix_',Fn, '_Run',run_idx, '.txt');

		if (!file.exists(factorFn)){
			#print(paste0(outputFn, ' does not exist'))
			next
		}
		FactorM     = read.table(factorFn, sep='\t');
		factor_corr = norm(cor(FactorM), 'F');

		F_all[[iii]] = FactorM;
		fcorr_all[[iii]] = factor_corr;
		iii = iii + 1;
		if(iii > 2){
			coph_test = c(coph_test, cophenet(F_all))
		}
	}


	if(length(F_all) > 1){
		coph = cophenet(F_all);
		cat('\n')
		print(paste(rep("#",30), collapse = ''));
		print(paste0(length(F_all), ' runs'));	
		print(coph_test);
		print(paste0('K = ', (K), '; alpha1 = ', (alpha1),'; lambda1 = ', (lambda1)));
		print(paste0('Coph = ', coph, "; mean fcorr = ", mean(unlist(fcorr_all))));
		print(paste(rep("#",30), collapse = ''))
		cat('\n');
		return(c(K, alpha1, lambda1, coph, mean(unlist(fcorr_all))))
	}else{
		print(paste0(Fn, ': not enough runs available'))
		return(NULL)
	}

}


result = NULL
for(K in c(30,35,40)){
	for(alpha1 in c(10,50, 100)){
		for(lambda1 in c(100,200,300,400,500)){
			rowi = collect_results(K, alpha1, lambda1)
			if(!is.null(rowi)){
				result = rbind(result, rowi)
			}
		}
	}
}


result = as.data.frame(result)
colnames(result) = c("K", "alpha1", "lambda1", "coph", "correlation")
result = result[order(result$coph), ]
write.table(result, '2_choose_hyperparameters_output.txt', sep='\t', quote = F, row.names = F)


