library(flashr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Rmosek)

library(ashr)
library(flashr)

suppressWarnings(library(readr))
source('sn_spMF/readIn.R')
source('sn_spMF/plot_factor_matrix.R')

#install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch",
#                 configure.vars = "PKG_MOSEKHOME=/Users/Yuan/Documents/mosek/8/tools/platform/osx64x86 PKG_MOSEKLIB=mosek64")
#mosek_attachbuilder(what_mosek_bindir='/Users/Yuan/Documents/mosek/8/tools/platform/osx64x86/bin', 
#                    pos=2L, name="Rmosek:builder", warn.conflicts=TRUE)
#install.rmosek()

xfn = 'data/test_data_X.txt'
wfn = 'data/test_data_W.txt'
rows_to_read = NULL
Data = readIn(xfn, wfn, rows_to_read);
X = Data[['X']];
W = Data[['W']];
Z_score = X * W ;


### From Gao: 
### https://gaow.github.io/mnm-gtex-v8/analysis/mashr_flashr_workflow.html#flashr-prior-covariances


my_init_fn <- function(Y, K = 1) {
  ret = flashr:::udv_si(Y, K)
  pos_sum = sum(ret$v[ret$v > 0])
  neg_sum = -sum(ret$v[ret$v < 0])
  if (neg_sum > pos_sum) {
    return(list(u = -ret$u, d = ret$d, v = -ret$v))
  } else
    return(ret)
}

## flashr       * 0.6-7     2020-01-10 [1] Github (stephenslab/flashr@ae44e7d)
## ashr           2.2-39    2020-01-10 [1] Github (stephens999/ashr@e8a7abc)
## ebnm         * 0.1-24    2020-01-06 [1] Github (stephenslab/ebnm@308cb8a) 
## mixsqp         0.2-2     2019-10-16 [1] CRAN (R 3.5.2)      

flash_pipeline = function(data, f_prior = '+uniform', Kmax = 50, 
                          outputdir = './output', saveFn = 'flashr', 
                          draw_factor = T, ...) {
  ## current state-of-the art
  ## suggested by Jason Willwerscheid
  ## cf: discussion section of
  ## https://willwerscheid.github.io/MASHvFLASH/MASHvFLASHnn2.html
  ebnm_fn = "ebnm_ash"
  ebnm_param = list(l = list(mixcompdist = "normal",
                             optmethod = "mixSQP"),
                    f = list(mixcompdist = f_prior,
                             optmethod = "mixSQP"))

  fl_g <- flashr:::flash_greedy_workhorse(data,
                                          var_type = "constant",
                                          ebnm_fn = ebnm_fn,
                                          Kmax = Kmax,
                                          ebnm_param = ebnm_param,
                                          init_fn = "my_init_fn",
                                          stopping_rule = "factors",
                                          tol = 1e-3,
                                          verbose_output = "odF")
  
    fl_b <- flashr:::flash_backfit_workhorse(data,
                                             f_init = fl_g,
                                             var_type = "constant",
                                             ebnm_fn = ebnm_fn,
                                             ebnm_param = ebnm_param,
                                             stopping_rule = "factors",
                                             tol = 1e-3,
                                             verbose_output = "odF")

    ## save the factor matrix
    flash_f = fl_b$ldf$f
    rownames(flash_f) = colnames(X)

    factorFn = paste0(outputdir, '/', saveFn, '.txt')
    write.table(flash_f, factorFn, sep='\t', row.names = T, quote = F)
    if(draw_factor){
        plot_factor_matrix(factorFn)
        plot_factor(factorFn)
    }
}




## flashr with default setting
fmodel = flash(as.matrix(Z_score), Kmax = 23)
# save the factor matrix
flash_f = fmodel$ldf$f
rownames(flash_f) = colnames(X)

outputdir = 'output/'
factorFn = paste0(outputdir, '/flashr.txt')
write.table(flash_f, factorFn, sep='\t', row.names = T, quote = F)
plot_factor_matrix(factorFn)
plot_factor(factorFn)


## non-negative factors
flash_pipeline(as.matrix(Z_score), Kmax = 23, saveFn = 'flashr_NN')


## factors with mixed signs
flash_pipeline(as.matrix(Z_score), f_prior = 'normal', Kmax = 23, saveFn = 'flashr_backfit')




