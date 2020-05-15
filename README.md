
# Learn factors using sn-spMF

## Getting Started

'''
`sn_spMF/`: Functions to run constraint matrix factorization

`run_MF.R`:  Main function to run weighted sn-spMF
            <br>
            > Objective: $\min_{F,L} ||(X - LF') \odot W||_F^2 + \alpha_1|L|_1 + \lambda_1|F|_1$, $F$ is non-negative
            <br>
            > Input: $K$: Rank of $F$, or number of factors;
                   $\alpha_1$: l1 penalty for $L$;
                   $\lambda_1$: l1 penalty for $F$.
		           (These parameters should be selected using `run_MF_train_coph.R`)
            <br>
            > Output: the factor matrix $F$
            <br>
            > Note: in R package "penalized", forced the penalization step to run with no standardization by setting weights = 1. This is specifically designed for eQTL effect sizes across tissues, since the strength of effect sizes provides information and should not be standardized.
'''


'''
`fit_F.R`: Update factors
            <br>
            > Objective: $\min_F ||(X - LF') \odot W||_F^2 + \lambda_1|F|_1$, $F$ is non-negative

`fit_L.R`: Update $L$
            <br>
            > Objective: $\min_L ||(X - LF') \odot W||_F^2 + \alpha_1|L|_1$
'''



##  Map eQTLs to factors using weighted linear regression

'''
`mapping/`: Functions to map eQTLs to factors

`run_linearReg.R`: Run weighted linear regression
'''
