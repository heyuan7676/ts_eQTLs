################################################################################################################################
## Update the factor matrix F during the alternative optimization
##
## Input:
##		- X: the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
##		- L: learned loading matrix  (N x K, N is the number of data points, K is the number of factors)
##		- W: the weight matrix, same size as in X
##		- option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
##
## Return:
##		- A non-negative factor matrix that minimize_F ||(X - LF') .* W||_F^2 + lambda1*|F|_1, F is non-negative
##
## Example of usage:
##
## source('../simulation/Generate_input.R');
## data = generate_input(tau = tau);
## L = data[[2]];
## X = data[[3]];
## W = data[[4]];
## option = list();
## option[['lambda1']] = 0.1;
## F_predict = fit_F(X, W, L, option);
##
################################################################################################################################

suppressWarnings(library(penalized))

fit_F <- function(X, W, L, option){
	tStart   = Sys.time();
	FactorM  = NULL;

	## fit each factor one by one -- because of the element-wise multiplication from weights!
	for(col in seq(1, ncol(X))){
		x = X[, col];
		w = W[, col];

		## weight the equations on both sides
		xp = w * x;
		Lp = w * L;

		## Fit: xp = L %*% f with l1 penalty on f -- |lambda1 * f|
		dat_i = as.data.frame(cbind(xp, Lp));
		colnames(dat_i) = c('X', paste0('F', seq(1, ncol(Lp))));

		# unpenalized = ~0 - suppress the intercept
		# positive = TRUE  - constrains the estimated regression coefficients of all penalized covariates to be non-negative
		# lambda2=1e-10    - avoid singularity in fitting
		fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
						unpenalized = ~0, lambda1 = option[['lambda1']], lambda2=1e-10,
						positive = T, standardize = F, trace = F)
		f = coef(fit, 'all')
		FactorM = rbind(FactorM, f);
	}

	tEnd = Sys.time();
	#print(paste0('Updating Factor matrix takes ', round((tEnd - tStart)/60, 2), 'min'));

	return(FactorM)

}
