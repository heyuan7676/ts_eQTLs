##############################################################################################################
## Compute objective: ||(X - LF') .* W||_F^2 + alpha1*|L|_1 + lambda1*|F|_1
##############################################################################################################

compute_obj <- function(X, W, L, FactorM, option){
	

	residual = W * (X - L %*% t(FactorM));
	Residual_penalty = sum(sum(residual^2));
	#Residual_penalty = sum(sum(residual^2)) / (2 * nrow(L));

	# the l2 penalty is added in the fitting step to avoid singularity, and is accounted for here
	L_penalty = option[['alpha1']]* sum(abs(L)) + 1e-10 * sum( L ^2);
	FactorM_penalty = option[['lambda1']] * sum(abs(FactorM)) + 1e-10 * sum(FactorM ^ 2);

	obj = Residual_penalty + FactorM_penalty + L_penalty;
	#print(paste0('obj decomposition: ', (Residual_penalty), '; ', (L_penalty), '; ', (FactorM_penalty)));

	return(obj)
}


