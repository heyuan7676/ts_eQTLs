script_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/Revision_Zscore/scripts/0_loadings/'
source(paste0(script_dir, "readin_YWX.R"))
source(paste0(script_dir, "run_linearReg.R"))
source(paste0(script_dir, "save_results.R"))

args <- commandArgs(TRUE)
FM_fn = args[1]
startidx = args[2]

## read in data
inputprefix = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1'
data  = readin_data(FM_fn, inputprefix, startidx = startidx)
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
save_results(Betas, pValues, FM_fn, startidx)



