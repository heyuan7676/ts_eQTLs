######################################
## Wrapper of updating L and F
## Input: X and W (both N x T)
## Output: L (N x K) and F (T x K)
######################################
 
source('compute_obj.R')
source('fit_L.R')
source('fit_F.R')


Update_FL <- function(X, W, option){
	# number of features - to avoid using T in R
	D      = dim(X)[2]
	tStart = Sys.time()

	#Random initialization
	FactorM   = matrix(runif(D*(option[['K']] - 1)), nrow = D);
	FactorM = cbind(rep(1, D), FactorM);

	# First round of optimization
	L = fit_L(X, W, FactorM, option);
	obj = c(compute_obj(X, W, L, FactorM, option));
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
		break
	}


	for (iii in seq(1, option[['iter']])){
		## update F
		obj_F_old  = obj[length(obj)];
        	FactorM = fit_F(X, W, L, option);

		# if number of factors decrease because of empty factor, the change in ||F||_F = 1
		if(dim(FactorM)[2] != dim(F_old)[2]){
			F_change = 1
		}else{
			F_change = norm(FactorM - F_old, 'F')
		}
		F_old = FactorM;
		obj_F  = compute_obj(X, W, L, FactorM, option);

		# if F is empty, stop
		non_empty_f = which(apply(FactorM, 2, function(x) sum(x!=0)) > 0)
        	if(length(non_empty_f) == 0){
                	print('Finished');
                	print('F is completely empty, lambda1 too large')
                	F_sparsity = 1;
                	L_sparsity = 1;
                	factor_corr = 1;
                	Nfactor = 0;
                	break
		}

		## update L
		obj_L_old = obj_F;
        	L  = fit_L(X, W, FactorM, option);
		obj_L = compute_obj(X, W, L, FactorM, option);

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
                	break
		}
		
		# align the two matrices	
		non_empty = intersect(non_empty_l,non_empty_f);
		L = L[, non_empty];
		FactorM  = FactorM[, non_empty];

		# collect sparsity in L and F
		L_sparsity = sum(abs(L) < 1e-5) / dim(L)[2] / dim(L)[1];
		F_sparsity = sum(abs(FactorM) < 1e-5) / dim(FactorM)[2] / dim(FactorM)[1];
		Nfactor = length(non_empty);

		# change in the objective -- should always be negative
		obj_change = (obj_F - obj_F_old) + (obj_L - obj_L_old);

		if(option[['disp']]){
			print(paste0('Iter', iii, ':'))
        		print(paste0('Loading Sparsity = ', L_sparsity, '; Factor sparsity = ', F_sparsity, '; ', Nfactor, ' factors remain'));
			print(paste0('Objective change for updating F = ', obj_F - obj_F_old, '; for updating L = ', obj_L - obj_L_old, '; Total change = ', obj_change));
			print(paste0('Factor matrix Frobenius norm change = ', F_change))
			cat('\n')
		}


		# converge if: 1). Change in ||F||_F is small. ie. The factor matrix is stable. 2). reached maximum number of iterations
		if(option[['toobj']] >= F_change){
                	print(paste0('Converged at itertaion ', (iii)));
                	break
		}
		if(iii == option[['iter']]){
        		print('Reached maximum iteration.');
        		break
		}

		old_change = obj_change;
		obj = c(obj, obj_L);
	}

	tEnd = Sys.time()
	print('Total time used for updating')
	print(tEnd - tStart);

	# return L, F, sparsity in L and F, number of factors -- could be different from K!
	return(list(FactorM, L, L_sparsity, F_sparsity, Nfactor))
}
