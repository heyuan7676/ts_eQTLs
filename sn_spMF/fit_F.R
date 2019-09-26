##############################################################################################################
## Update factors
## Objective: minimize_F ||(X - LF') .* W||_F^2 + lambda1*|F|_1, F is non-negative
##############################################################################################################

#library(penalized)

fit_F <- function(X, W, L, option){
        tStart   = Sys.time();
        FactorM  = NULL;

	## fit each factor one by one -- because of the element-wise multiplication from weights!
        for(col in seq(1, dim(X)[2])){
                x = X[, col];
                w = W[, col];

		## weight the equations on both sides
		xp = w * x;
        	Lp = w * L;

                ## Fit: xp = L %*% f with l1 penalty on f -- |lambda1 * f|
                dat_i = as.data.frame(cbind(xp, Lp));
                colnames(dat_i) = c('X', paste0('F', seq(1, dim(Lp)[2])));

		# unpenalized = ~0 - suppress the intercept
		# positive = TRUE  - constrains the estimated regression coefficients of all penalized covariates to be non-negative
		# lambda2=1e-10    - avoid singularity in fitting
                fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, dim(Lp)[2]))], data=dat_i,
                                unpenalized = ~0, lambda1 = option[['lambda1']], lambda2=1e-10,
                                positive = T, standardize = F, trace = F)
		f = coef(fit, 'all')
                FactorM = rbind(FactorM, f);
        }

        tEnd = Sys.time();
        #print('Updating Factor matrix: ');
        #print(tEnd - tStart)

	colnames(FactorM) = seq(1, dim(FactorM)[2])
        return(FactorM)

}
