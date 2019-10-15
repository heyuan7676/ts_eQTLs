**Code for learning factors using sn-spMF and mapping eQTLs to factors using weighted LR**

In sn_spMF/: Functions to run constraint matrix factorization

`run_MF.R`: Main function to run weighted sn-spMF \
            Objective: minimize_{F,L} ||(X - LF') .* W||_F^2 + alpha1*|L|_1 + lambda1*|F|_1, F is non-negative \
            Input: K: Rank of F (number of factors), alpha1: l1 penalty for L, lambda1: l1 penalty for F. \
		           (These parameters should be selected using run_MF_train_coph.R) \
            Output: the factor matrix F (FactorM) \
            Note: in R package "penalized", forced the penalization step to run with no standardization by setting weights = 1. This is specifically designed for eQTL effect sizes across tissues, since the strength of effect sizes provides information and should not be standardized. \

`fit_F.R`: Update factors\
            Objective: minimize_F ||(X - LF') .* W||_F^2 + lambda1*|F|_1, F is non-negative

`fit_L.R`: Update L \
            Objective: minimize_L ||(X - LF') .* W||_F^2 + alpha1*|L|_1


sn_spMF/run_MF_train_coph.R: The wrapper for parameter tuning in sn_spMF\
    -- readIn.R: Read in data. Please modify to fit\
    -- Update_FL.R: Implement alternating least squares (ALS) with gradient descent\
        -- fit_F.R: Update factors\
        -- fit_L.R: Update loadings\
        -- compute_obj.R: Compute the objective\
    -- compute_cophenet.R: Compute cophenetic coefficient after multiple runs.\

sn_spMF/run_MF.R: The wrapper for running sn-spMF.\ 
    -- readIn.R: Read in data. Please modify to fit\
    -- Update_FL.R: Implement alternating least squares (ALS) with gradient descent\
	-- fit_F.R: Update factors\
	-- fit_L.R: Update loadings\
	-- compute_obj.R: Compute the objective \


mapping/run_linearReg.R - mapping eQTLs to the factors
