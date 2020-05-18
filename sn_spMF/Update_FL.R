######################################
## Wrapper of updating L and F
## Input: X and W (both N x T)
## Output: L (N x K) and F (T x K)
######################################
 
source('sn_spMF/compute_obj.R')
source('sn_spMF/fit_L.R')
source('sn_spMF/fit_F.R')


Update_FL <- function(X, W, option){
	# number of features - to avoid using T in R
	D      = ncol(X)
	tStart0 = Sys.time()

	#Random initialization
	print('Random initialization...')
	FactorM   = matrix(runif(D*(option[['K']] - 1)), nrow = D);
	FactorM = cbind(rep(1, D), FactorM);

	# First round of optimization
	print('Start optimization ...')
	L = fit_L(X, W, FactorM, option);
	objective = c(NA, compute_obj(X, W, L, FactorM, option));
	objective_change = c(1, 1);
	old_change = 1;
	F_old = FactorM; 

	# if L is empty, stop
	non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
	if(length(non_empty_l) == 0){
		print('Finished');
		print('L is completely empty, alpha1 too large')
		FactorM = NULL;
		F_sparsity = 1;
		L_sparsity = 1;
		factor_corr = 1;
		Nfactor = 0;
		return()
	}


	for (iii in seq(1, option[['iter']])){

		## update F
		FactorM = fit_F(X, W, L, option);

		# if number of factors decrease because of empty factor, the change in ||F||_F = 100
		if(ncol(FactorM) != ncol(F_old)){
			F_change = 100
		}else{
			F_change = norm(FactorM - F_old, 'F') / ncol(FactorM)
		}
		F_old = FactorM;

                non_empty_f = which(apply(FactorM, 2, function(x) sum(x!=0)) > 0)
                if(length(non_empty_f) == 0){
                        print('Finished');
                        print('F is completely empty, lambda1 too large')
                        F_sparsity = 1;
                        L_sparsity = 1;
                        factor_corr = 1;
                        Nfactor = 0;
                        return();
                }
		colnames(FactorM) = seq(1, ncol(FactorM));

		## update L
        	L  = fit_L(X, W, FactorM, option);

		# if L is empty, stop
        	non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
        	if(length(non_empty_l) == 0){
                        print('Finished');
                        print('L is completely empty, alpha1 too large')
                	FactorM = NULL;
                	F_sparsity = 1;
                	L_sparsity = 1;
                	factor_corr = 1;
                	Nfactor = 0;
                	return()
		}
		colnames(L) = seq(1, ncol(L));

		# align the two matrices	
		non_empty = intersect(non_empty_l,non_empty_f);
		L = as.matrix(as.data.frame(L[, non_empty]));
		FactorM  = as.matrix(as.data.frame(FactorM[, non_empty]));


		# collect sparsity in L and F
		L_sparsity = sum(abs(L) < 1e-5) / ncol(L) / nrow(L);
		F_sparsity = sum(abs(FactorM) < 1e-5) / ncol(FactorM) / nrow(FactorM);
		Nfactor = length(non_empty);

		# change in the objective -- should always be negative
		obj_updated = compute_obj(X, W, L, FactorM, option);
		obj_change = obj_updated - objective[length(objective)];
		objective = c(objective, obj_updated);
		objective_change = c(objective_change, obj_change);

		if(option[['disp']]){
			cat('\n')
			print(paste0('Iter', iii, ':'))
			print(paste0('Objective change = ', obj_change))
			print(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', F_change));
			print(paste0('Loading Sparsity = ', L_sparsity, '; Factor sparsity = ', F_sparsity, '; ', Nfactor, ' factors remain')); 
			cat('\n')
		}


		# converge if: 1). Change of the values in factor matrix is small. ie. The factor matrix is stable. 2). Change in the objective function becomes small; 3). reached maximum number of iterations
		if(option[['convF']] >= F_change){
                	print(paste0('Factor matrix converges at itertaion ', iii));
                	break
		}

		oc1 = abs(objective_change[length(objective_change)])
		if(option[['convO']] >= oc1){
                        print(paste0('Objective function converges at itertaion ', iii));
                        break
                }
		if(iii == option[['iter']]){
        		print('Reached maximum iteration.');
        		break
		}

	}

	tEnd0 = Sys.time()
	cat('\n')
	print('Total time used for optimization: ');
	print(tEnd0 - tStart0);

	# return L, F, sparsity in L and F, number of factors -- could be different from K!
	return(list(FactorM, L, L_sparsity, F_sparsity, Nfactor, objective))
}
