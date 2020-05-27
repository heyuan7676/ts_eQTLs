source("mapping/readin_data.R")
source("mapping/run_linearReg.R")
source("mapping/adj_pvalue.R")
source("mapping/Remove_compensate_pairs.R")

fit_significant_hits <- function(F_C, FM_fn_save, xfilename = 'simulation/input/Input_tau1000_seed1_X.txt', 
				wfilename = 'simulation/input/Input_tau1000_seed1_W.txt', 
				mappingDir = 'simulation/output/mapping/'){

	factorDir = 'notused'
	FM_fn_tp = 'notused'

	## read in data
	data  = readin_YWX(FM_fn_tp, xfn=xfilename, wfn=wfilename, factor_output_dir = factorDir, F_C = F_C)
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
	adj_pvalue(Betas, pValues, FM_fn_save, mappingDir)
	Data = readin_XWBF(FM_fn_save, xfn=xfilename, wfn=wfilename, 
		   mapping_output_dir = mappingDir, factor_output_dir = factorDir, F_C = F_C)
	X = Data[[1]]
	W = Data[[2]]
	B = Data[[3]]
	FactorM = Data[[4]]


	return(as.matrix(B))
}
