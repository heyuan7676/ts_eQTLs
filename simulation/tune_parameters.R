source("sn_spMF/cophenet.R")
source("sn_spMF/collect_results.R")
source("sn_spMF/fit_L.R")
source("sn_spMF/readIn.R")

suppressWarnings(library(optparse))

option_list = list(make_option(c("-O", "--outputdir"), type = "character", default='simulation/output/', help="output directory", metavar="character"))
opt = parse_args(OptionParser(option_list=option_list))


result = NULL
for(K in seq(5,9)){
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

bnm = basename(opt$outputdir);
write.table(result, paste0('simulation/output/', bnm, '_choose_paras.txt'), sep='\t', quote = F, row.names = F)

