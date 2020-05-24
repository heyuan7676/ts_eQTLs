suppressWarnings(library(plyr))


library(optparse) 
source("sn_spMF/sn_spMF.R")

option_list = list(make_option(c("-O", "--outputdir"), type = "character", default='output/', help="output directory", metavar="character"),
				   make_option(c("-k", "--k"), type="integer", default=30,help="Number of factors"),
				   make_option(c("-a", "--alpha1"), default = 0.01, type="double", help="l1 penalty for L"),
				   make_option(c("-l", "--lambda1"), default = 0.1, type="double", help="l1 penalty for F"))


find_optimal <- function(outputdir, K, alpha1, lambda1){

	# read in the arguments
	Fn = paste0('K', K, '_a1', alpha1, '_l1', lambda1); 
	outputdir = file.path(outputdir, paste0('sn_spMF_', Fn));

	min_obj = Inf
	for(run_idx in seq(1,30)){
		outputFn = paste0(outputdir, '/sn_spMF_', Fn, '_Run',run_idx, '.RData');
		if (!file.exists(outputFn)){
			print(paste0(outputFn, ' does not exist'))
			next
		}

		load(outputFn);
		if(objective[length(objective)] < min_obj){ 
			fn_optimal = outputFn
			min_obj = objective[length(objective)]
		}
	}
	
	print(paste0('Optimal solution: ', gsub(".RData", "", fn_optimal)))
}


opt = parse_args(OptionParser(option_list=option_list))
find_optimal(opt$outputdir, opt$k, opt$alpha1, opt$lambda1)



