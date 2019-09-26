source("readin_YWX.R")
source("run_linearReg.R")
source("save_results.R")

args <- commandArgs(TRUE)
FM_fn    = args[1]
inputprefix   = args[2]
LDr2 = args[3]
startidx = args[4]

## read in data
data  = readin_data(FM_fn, paste0(inputprefix, '_', LDr2), startidx)
Genes = data[[1]]
SNPs  = data[[2]]
dataPoints = data[[3]]
comp = data[[4]]

## run regression
result  = run_regression(dataPoints, comp)
Betas   = result[[1]]
pValues = result[[2]]

rownames(pValues) = paste(Genes, SNPs)
rownames(Betas)   = paste(Genes, SNPs)

## save results
save_results(Betas, pValues, FM_fn, paste0('_LD', LDr2), startidx)



