suppressWarnings(library(ssvd))
suppressWarnings(library(softImpute))
suppressWarnings(library(PMA))
suppressWarnings(library(readr))

source('sn_spMF/readIn.R')
source('sn_spMF/plot_factor_matrix.R')
source('simulation//utils.R')

xfn = 'data/test_data_X.txt'
wfn = 'data/test_data_W.txt'
rows_to_read = NULL
Data = readIn(xfn, wfn, rows_to_read);
X = Data[['X']];
W = Data[['W']];
Z_score = X * W ;

rankK = 23

###################################################
### softImpute
###################################################

l1 = 50
while(TRUE){
	softimputes = softImpute(as.matrix(Z_score),lambda = 2 * l1, rank.max = rankK)
	softimpute_f = softimputes$v
        if (dim(softimpute_f)[2] < 23){
                break
        }else{
		l1 = 2*l1
	}
}


while(TRUE){
        softimputes = softImpute(as.matrix(Z_score),lambda = l1 + 10, rank.max = rankK)
        softimpute_f = softimputes$v
        if (dim(softimpute_f)[2] < 23){
                break
        }else{
                l1 = l1 + 10
        }
}

softimputes = softImpute(as.matrix(Z_score),lambda = l1, rank.max = rankK)
softimpute_f = softimputes$v
softimpute_f = as.data.frame(softimpute_f)
rownames(softimpute_f) = colnames(X)

# save
outputdir = 'output/'
factorFn = paste0(outputdir, '/softeImput.txt')
write.table(softimpute_f, factorFn, sep='\t', row.names = T, quote = F)
plot_factor_matrix(factorFn)
plot_factor(factorFn)
