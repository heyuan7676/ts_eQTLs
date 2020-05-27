suppressWarnings(library(plyr))
suppressWarnings(library(optparse))

option_list = list(make_option(c("-O", "--outputdir"), type = "character", default='output/', help="output directory", metavar="character"),
		   make_option(c("-f", "--savefn"), type = "character", default='choose_para_preliminary.txt', help="filename to save the output", metavar="character"))
opt = parse_args(OptionParser(option_list=option_list))


collect_results_preliminary <- function(outputdir, K, alpha1, lambda1){
        # read in the arguments
        Fn = paste0('K', K, '_a1', alpha1, '_l1', lambda1);
        outputdir = file.path(outputdir, paste0('sn_spMF_', Fn));
	outputFn = paste0(outputdir, '/sn_spMF_', Fn, '.RData');
	load(outputFn);

	return(c(K, alpha1, lambda1, ncol(FactorM), L_sparsity, F_sparsity))
                                                                       
}



result = NULL
for(K in c(10, 15, 20)){
        for(alpha1 in c(1,10,100,500)){
                for(lambda1 in c(1,10,100,500)){
                        rowi = collect_results_preliminary(opt$outputdir, K, alpha1, lambda1)
                        if(!is.null(rowi)){
                                result = rbind(result, rowi)
                        }
                }
        }
}

result = as.data.frame(result)
colnames(result) = c("K", "alpha1", "lambda1", "nfactor", "L_sparsity", "F_sparsity")
write.table(result, paste0(opt$outputdir, opt$savefn), sep='\t', quote = F, row.names = F)


