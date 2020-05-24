source("mapping/readin_data.R")
source("mapping/run_linearReg.R")
source("mapping/adj_pvalue.R")
source("mapping/Remove_compensate_pairs.R")

library(optparse)

option_list = list(
        make_option(c("-x", "--xfilename"), type = "character", default='data/test_data_X.txt', help="input filename for X", metavar="character"),
        make_option(c("-w", "--wfilename"), type = "character", default='data/test_data_SE.txt', help="input filename for W", metavar="character"),
        make_option(c("-f", "--factorMname"), type = "character", default='sn_spMF_K15_a1100_l130', help="factor matrix filename", metavar="character"),
	make_option(c("-d", "--factorDir"), type = "character", default='output/', help="output directory for the factor matrix", metavar="character"),
	make_option(c("-m", "--mappingDir"), type = "character", default='output/mapping/', help="output directory for the loadings", metavar="character"))

opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$mappingDir, showWarnings = F)


FM_fn = opt$factorMname

## read in data
data  = readin_YWX(FM_fn, xfn=opt$xfilename, wfn=opt$wfilename, factor_output_dir = opt$factorDir)
Genes = data[[1]]
SNPs  = data[[2]]
dataPoints = data[[3]]
compositions = data[[4]]
tissues = data[[5]]

## run regression
result  = run_regression(dataPoints, tissues, compositions)
Betas   = result[[1]]
pValues = result[[2]]

rownames(pValues) = paste(Genes, SNPs)
rownames(Betas)   = paste(Genes, SNPs)

## Perform BH-correction
adj_pvalue(Betas, pValues, FM_fn, opt$mappingDir)


## remove compensave effects
Data = readin_XWBF(FM_fn, xfn=opt$xfilename, wfn=opt$wfilename, 
		   mapping_output_dir = opt$mappingDir, factor_output_dir = opt$factorDir)
X = Data[[1]]
W = Data[[2]]
B = Data[[3]]
FactorM = Data[[4]]

B_correct = correct_B(X, W, B, FactorM)

## save results
outFN = paste0('output/mapping/mapping_', FM_fn, '_Loadings_beta_alpha0.05_corrected.txt')
print(outFN)
write.table(B_correct, outFN, sep='\t', quote = F)


print('Proportion of eQTLs that load on each factor')
print(apply(B_correct, 2, function(x) sum(x!=0)) / nrow(B_correct))
