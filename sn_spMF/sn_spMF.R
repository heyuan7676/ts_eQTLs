##############################################################################################################################
## Main function to run weighted sn-spMF 
##
## Objective: minimize_{F,L} ||(X - LF') .* W||_F^2 + alpha1*|L|_1 + lambda1*|F|_1, F is non-negative
##
## Input: 
##	- xfn: filename of the input matrix for decomposition
##	- wfn: filename of the weight matrix, can be NULL if there is no weight matrix
## 	- K:   Rank of F (or number of factors) for initialization
##	- alpha1: l1 penalty for loading matrix (L)
##	- lambda1: l1 penalty for factor matrix (F)
## 	- run_idx: indicate whether perform one run in one process, or perform 30 runs in one process (run_idx = -1 by default).
##			If the dataset has moderate size (like 1k data points), it is reasonable to run the 30 runs in one script to control the number of scripts running
##			Otherwise if the dataset is large (ie. 10k data points), you may want to run each run in parallel to speed up
##	- ... (other arguments)
##
## The function does:
##	- save the learned factor matrix and corresponding arguments
##	- plot the learned factor matrix (if draw_factor is TRUE)
##
## Example of usage: see run_MF.R
## 
##############################################################################################################################

source("sn_spMF/readIn.R")
source("sn_spMF/Update_FL.R")
source("sn_spMF/cophenet.R")
source("sn_spMF/plot_factor_matrix.R")



sn_spMF <- function(xfn, wfn, K, alpha1, lambda1, outputdir, 
			iterations = 5, converged_F_change = 0.1, converged_obj_change = 0.1,
			verbose = 1, compute_cophenet = 0, run_idx = -1, 
			draw_factor = 1, rows_to_read = NULL){

	dir.create(outputdir, showWarnings = FALSE)

	# read in the data points
	## X - data matrix to decompose. N x T, where N is the number of data points, T is the number of features
	## W - matrix of weigths for each data point. N x T
	Data = readIn(xfn, wfn, rows_to_read);
	X = Data[['X']];
	W = Data[['W']];

	if(K > ncol(X)){
		print('ERROR: Number of factors should not exist number of features')
		quit()
	}
	print(paste0('Fit sn-spMF for data of ', nrow(X), ' data points and ', ncol(X), ' features'));

	# read in the arguments
	option = list()
	option[['alpha1']]  = alpha1;
	option[['lambda1']] = lambda1;
	option[['K']]       = K;
	option[['iter']]  = iterations;
	option[['convF']] = converged_F_change;
	option[['convO']] = converged_obj_change;
	option[['disp']]  = verbose == 1;

	Fn = paste0('K', K, '_a1', alpha1, '_l1', lambda1); 
        outputdir = file.path(outputdir, paste0('sn_spMF_', Fn))
        dir.create(outputdir,  showWarnings = FALSE)

	# Run the model multiple times with random initialization - to compute the stability of the decomposed matrices
	if(compute_cophenet == 1){
		print('Tune the hyper-parameters ... ')

		## If the dataset has moderate size (like 1k data points), it is reasonable to run the 30 runs in one script
		## Otherwise if the dataset is large (ie. 10k data points), you may want to run each run in parallel to speed up
		if(run_idx < 0){
			runs = seq(1,30)
		}else{
			runs = c(run_idx)
		}

		for(run_idx in runs){
			cat('\n')
			print(paste(rep("#",30), collapse = ''))
			print(paste0('## Run', run_idx));
			print(paste(rep("#",30), collapse = ''))
			cat('\n')
			run_one_iter(X, W, option, outputdir, paste0(Fn, '_Run',run_idx), draw_factor);
		}

	}else{
		run_one_iter(X, W, option, outputdir, Fn, draw_factor);
	}

}




run_one_iter <- function(X, W, option, outputdir, outputFn, draw_factor){
        Run_iter = Update_FL(X, W, option);

        if(is.null(Run_iter)){
                return()
        }


        FactorM = Run_iter[[1]];
        factor_corr = norm(cor(FactorM), 'F');
	LoadingM = Run_iter[[2]];
        L_sparsity = Run_iter[[3]];
        F_sparsity = Run_iter[[4]];
        objective = Run_iter[[6]];

	alpha1 = option[['alpha1']]
	lambda1 = option[['lambda1']]
	K = option[['K']]

        cat('\n')
        print('Optimization finished')
        print(paste0('K = ', K, '; alpha1 = ', alpha1,'; lambda1 = ', lambda1));
        print(paste0('Sparsity in Factor matrix = ', (F_sparsity),'; Sparsity in L = ', (L_sparsity), '; '))
        print(paste0((Run_iter[[5]]), ' factors remain; correlation between factors = ', (factor_corr)));


        # store results
	rownames(FactorM) = colnames(X);
        Fn = paste0(outputdir, '/sn_spMF_',outputFn, '.RData');
        save(alpha1, lambda1, K, FactorM, objective, file = Fn);

        factorFn = paste0(outputdir, '/sn_spMF_FactorMatrix_',outputFn, '.txt');
        write.table(FactorM, factorFn, sep='\t', quote = F);

        if(draw_factor){
                plot_factor_matrix(factorFn)
                plot_factor(factorFn)
        }
}
