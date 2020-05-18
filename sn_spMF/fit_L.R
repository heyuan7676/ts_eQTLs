##############################################################################################################
## Update L
## Objective: minimize_L ||(X - LF') .* W||_F^2 + alpha1*|L|_1
## Similar to fit_F.R
##############################################################################################################

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


