##############################################################################################################
## Update L
## Objective: minimize_L ||(X - LF') .* W||_F^2 + alpha1*|L|_1
## Similar to fit_F.R
##############################################################################################################

#library(penalized)

one_fit <- function(x, w, FactorM, option){
                xp = w * x;
                FactorMp = diag(w) %*% FactorM;

                ## Fit: xp' = FactorMp %*% l with l1 penalty on l -- |alpha1 * l|
                dat_i = as.data.frame(cbind(t(xp), FactorMp));
                colnames(dat_i) = c('X', paste0('F', seq(1, dim(FactorMp)[2])));

                fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,dim(FactorMp)[2]))], data=dat_i,
                                unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                                positive = F, standardize = F, trace = F)
		l = coef(fit, 'all')
                return(l)
}


fit_L <- function(X, W, FactorM, option){

	L = NULL
        tS = Sys.time()
        for(row in seq(1,dim(X)[1])){
                x = X[row, ];
                w = W[row, ];
                l = one_fit(x, w, FactorM, option)
                L = rbind(L, l)
        }
        tE = Sys.time()
	#print('Updating Loading matrix: ');
        #print(tE - tS)

	colnames(L) = seq(1, dim(L)[2])
	return(L)
}
