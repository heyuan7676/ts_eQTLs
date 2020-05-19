################################################################################################################################
## Update the factor matrix L during the alternative optimization
##
## Input:
##              - X: the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
##              - F: learned factor matrix  (T x K, T is the number of features, K is the number of factors)
##              - W: the weight matrix, same size as in X
##              - option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
##
## Return:
##              - A loading matrix with mixed signs that minimize_F ||(X - LF') .* W||_F^2 + lambda1*|L|_1
##
## Example of usage:
##
## source('../simulation/Generate_input.R');
## data = generate_input(tau = tau);
## F = data[[1]];
## X = data[[3]];
## W = data[[4]];
## option = list();
## option[['lambda1']] = 0.1;
## L_predict = fit_L(X, W, F, option);
##
################################################################################################################################


fit_L <- function(X, W, FactorM, option){
	L = NULL
	tS = Sys.time()
	for(row in seq(1,nrow(X))){
		x = X[row, ];
		w = W[row, ];
		l = one_fit(x, w, FactorM, option);
		L = rbind(L, l);
	}
	tE = Sys.time();
	#print(paste0('Updating Loading matrix takes ', round((tE - tS)/60, 2), 'min'));

	return(L)
}



one_fit <- function(x, w, FactorM, option){
	xp = w * x;
	FactorMp = diag(w) %*% FactorM;

	# Fit: xp' = FactorMp %*% l with l1 penalty on l -- |alpha1 * l|
	dat_i = as.data.frame(cbind(t(xp), FactorMp));
	colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));

	fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
					unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
					positive = F, standardize = F, trace = F);
	l = coef(fit, 'all');
	return(l)
}


