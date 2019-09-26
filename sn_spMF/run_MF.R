##############################################################################################################################
## Main function to run weighted sn-spMF 
##
## Objective: minimize_{F,L} ||(X - LF') .* W||_F^2 + alpha1*|L|_1 + lambda1*|F|_1, F is non-negative
##
## Input: K - Rank of F / number of factors, alpha1- l1 penalty for L, lambda1 - l1 penalty for F.
## These parameters should be selected using run_MF_train_coph.R
##
## Output: the factor matrix F (FactorM)
##
## Note: in R package "penalized", forced the penalization step to run with no standardization by setting weights = 1.
##	 This is specifically designed for eQTL effect sizes across tissues, since the strength of effect sizes
##	 provides information and should not be standardized. 
##
##############################################################################################################################

#install.packages('optparse')
#install.packages('penalized')
#install.packages('plyr')

source("readIn.R")
source("Update_FL.R")
library(penalized)
library(optparse)

option_list = list(
	make_option(c("-i", "--inputdir"), type = "character", default='/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/', help="input directory", metavar="character"),
	make_option(c("-x", "--xfilename"), type = "character", default='test_data_X.txt', help="input filename for X", metavar="character"),
	make_option(c("-w", "--wfilename"), type = "character", default='test_data_W.txt', help="input filename for W", metavar="character"),
	make_option(c("-o", "--outputdir"), type = "character", default='/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/R_output/', help="output directory", metavar="character"),
	make_option(c("-k", "--k"), type="integer", default=20,help="Number of factors"),
	make_option(c("-a", "--alpha1"), type="double", help="l1 penalty for L"),
	make_option(c("-l", "--lambda1"), type="double", help="l1 penalty for F"),
	make_option(c("-t", "--iterations"), type="integer", default=5,help="Maximum number of iterations"))

opt = parse_args(OptionParser(option_list=option_list))

K=opt$k
alpha1=opt$alpha1
lambda1=opt$lambda1

inputdir=opt$inputdir
xfn=paste0(inputdir,opt$xfilename)
wfn=paste0(inputdir,opt$wfilename)
outputdir=opt$outputdir

## read in the data points
## X - data matrix to decompose. N x T, where N is the number of data points, T is the number of features
## W - matrix of weigths for each data point. N x T
## option - include K, alpha1, lambda 
Data = readIn(K, alpha1, lambda1, xfn, wfn);
X = Data[['X']];
W = Data[['W']];
option = Data[['option']];
option[['iter']]  = opt$iterations;
Fn_basename = Data[['Fn_basename']];

## run MF to learn the factors  
print(paste0('K=', (K), '; alpha1=', (alpha1),'; lambda1=', (lambda1)));
Run_iter = Update_FL(X, W, option);
FactorM = Run_iter[[1]]
factor_corr = norm(cor(FactorM), 'F')
L_sparsity = Run_iter[[3]]
F_sparsity = Run_iter[[4]]

print(paste0('Sparsity in Factor matrix =', (F_sparsity),'; Sparsity in L =', (L_sparsity), '; '))
print(paste0((Run_iter[[5]]), ' factors remain; ; correlation between factors = ', (factor_corr)));

# store results
Fn = paste0(outputdir, 'SparseMF_results_',Fn_basename, '.RData');
save(alpha1, lambda1, K, FactorM, file = Fn);

