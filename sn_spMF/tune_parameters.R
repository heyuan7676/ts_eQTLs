source("sn_spMF/cophenet.R")
source("sn_spMF/collect_results.R")
suppressWarnings(library(plyr))
suppressWarnings(library(optparse))

option_list = list(make_option(c("-O", "--outputdir"), type = "character", default='output/', help="output directory", metavar="character"),
		   make_option(c("-f", "--savefn"), type = "character", default='choose_para.txt', help="filename to save the output", metavar="character"))
opt = parse_args(OptionParser(option_list=option_list))


result = NULL
for(K in seq(10,20)){
	for(alpha1 in seq(1,10)){
		for(lambda1 in seq(1,10)){
			rowi = collect_results(opt$outputdir, K, 10*alpha1, 10*lambda1)
			if(!is.null(rowi)){
				result = rbind(result, rowi)
			}
		}
	}
}


result = as.data.frame(result)
colnames(result) = c("K", "alpha1", "lambda1", "coph", "correlation", "nfactor", "optimal_run", "optimal_obj")
result = result[order(result$coph), ]
write.table(result, paste0(opt$outputdir, opt$savefn), sep='\t', quote = F, row.names = F)


