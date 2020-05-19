##############################################################################################################################
## Wrapper to run sn-spMF 
##
## Example of usage: see run_MF.R
## Rscript run_MF.R -k 20 -a 10 -l 200 -t 10 -c 1 -r 1
##
##############################################################################################################################


#install.packages('optparse', repo="http://cran.rstudio.com/")
#install.packages('penalized', repo="http://cran.rstudio.com/")
#install.packages('plyr', repo="http://cran.rstudio.com/")

library(optparse) 
source("sn_spMF.R")

option_list = list(
	make_option(c("-x", "--xfilename"), type = "character", default='../data/test_data_X.txt', help="input filename for X", metavar="character"),
	make_option(c("-w", "--wfilename"), type = "character", default='../data/test_data_W.txt', help="input filename for W", metavar="character"),
	make_option(c("-O", "--outputdir"), type = "character", default='../output/', help="output directory", metavar="character"),
	make_option(c("-k", "--k"), type="integer", default=30,help="Number of factors"),
	make_option(c("-a", "--alpha1"), default = 0.01, type="double", help="l1 penalty for L"),
	make_option(c("-l", "--lambda1"), default = 0.1, type="double", help="l1 penalty for F"),
	make_option(c("-t", "--iterations"), type="integer", default=50,help="Maximum number of iterations"),
	make_option(c("-f", "--converged_F_change"), type="double", default=0.02,help="Change in factor matrix to call convergence"),
	make_option(c("-o", "--converged_obj_change"), type="double", default=1,help="Relative change in the objective function to call convergence"),
	make_option(c("-c", "--cophenet"), type="integer", default=0,help="Run the model multiple times to compute stability or not"),
	make_option(c("-r", "--runidx"), type="integer", default=3,help="The ith runs to compute cophenetic coefficient"),
	make_option(c("-v", "--verbose"), type="integer", default=1,help="Display the progress"))


opt = parse_args(OptionParser(option_list=option_list))
sn_spMF(opt$xfilename, 
		opt$wfilename, 
		opt$k, 
		opt$alpha1, 
		opt$lambda1, 
		outputdir = opt$outputdir, 
		iterations = opt$iterations, 
		converged_F_change = opt$converged_F_change,
		converged_obj_change = opt$converged_obj_change,
		verbose = opt$verbose,
		compute_cophenet = opt$cophenet, 
		run_idx = opt$runidx, 
		rows_to_read = NULL)
