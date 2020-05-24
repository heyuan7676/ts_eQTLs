source("sn_spMF/cophenet.R")
source("sn_spMF/collect_results.R")
suppressWarnings(library(plyr))
suppressWarnings(library(optparse))

option_list = list(make_option(c("-O", "--outputdir"), type = "character", default='output/', help="output directory", metavar="character"),
		   make_option(c("-f", "--savefn"), type = "character", default='choose_para_preliminary.txt', help="filename to save the output", metavar="character"))
opt = parse_args(OptionParser(option_list=option_list))



result = NULL
for(K in c(10, 15, 20)){
        for(alpha1 in c(1,10,100,500)){
                for(lambda1 in c(1,10,100,500)){
                        rowi = collect_results(opt$outputdir, K, alpha1, lambda1)
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


