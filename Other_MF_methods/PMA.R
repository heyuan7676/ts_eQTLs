suppressWarnings(library(ssvd))
suppressWarnings(library(softImpute))
suppressWarnings(library(PMA))
suppressWarnings(library(readr))

source('sn_spMF/readIn.R')
source('sn_spMF/plot_factor_matrix.R')
source('sn_spMF/utils.R')

xfn = 'data/test_data_X.txt'
wfn = 'data/test_data_W.txt'
rows_to_read = NULL
Data = readIn(xfn, wfn, rows_to_read);
X = Data[['X']];
W = Data[['W']];
Z_score = X * W ;

rankK = 23


###################################################
### PMA_cv1
###################################################

x = as.matrix(Z_score)
pma_f = NULL
for(k in seq(1,23)){
	cv.out <- PMD.cv(x, type="standard", sumabss=seq(0.2, 1, len=20))
	pmas = PMD(as.matrix(x), sumabs = cv.out$bestsumabs, sumabsu=NULL, sumabsv=NULL) ## sumabs: sparsity of v, sumabsu: sparsity of u
	x = x - pmas$u %*% pmas$d %*% t(pmas$v)
	pma_f = cbind(pma_f, pmas$v)
} 

# save
outputdir = 'output/'
factorFn = paste0(outputdir, '/PMA_cv1.txt')
write.table(pma_f, factorFn, sep='\t', row.names = T, quote = F)
plot_factor_matrix(factorFn)
plot_factor(factorFn)


###################################################
### PMA_cv2
###################################################

result = NULL
for(sv in seq(1,sqrt(49))){
	for(su in c(1,10,50, 100,120)){
		cv = tune_cor_PMA(Z_score, su, sv)
		result = rbind(result, c(sv, su, cv))
	}
}

ops = result[which.min(result[,3]),]
sv = ops[1]
su = ops[2]
pmas = PMD(as.matrix(Z_score), K=rankK, sumabs = NULL, sumabsu = su, sumabsv = sv)

# save
outputdir = 'output/'
factorFn = paste0(outputdir, '/PMA_cv2.txt')
write.table(pma_f, factorFn, sep='\t', row.names = T, quote = F)
plot_factor_matrix(factorFn)
plot_factor(factorFn)


