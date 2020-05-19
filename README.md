
# Learn latent factors using sn-spMF
We develop a constrained matrixfactorization model to learn patterns of tissue-sharing and tissue-specificity of eQTLs across 49 human tissuesfrom the Genotype-Tissue Expression (GTEx) project. The learned factors include patterns reflecting tissueswith known biological similarity or shared cell types, in addition to a dense factor representing a ubiquitousgenetic effect across all tissues

### Prerequisites
You need to install R/3.5.1 to run the scripts. R packages needed are:
```
install.packages('penalized')
```

## Running the tests
```
Rscript run_MF.R -k 10 -a 10 -l 100 -t 100
```

### Input files

```data/test_data_X.txt```: each row contains the effect size of an eQTL across tissues; the first two columns are gene names and SNP names for the eQTLs, and following columns are the features to learn patterns about, (tissues in the demo, can be time points in time-series data, or cells in single cell data). Missing data are presented as NA. Columns are seperated by '\t'. 

```
Gene	SNP	Adipose_Subcutaneous	Adipose_Visceral_Omentum	Adrenal_Gland	Artery_Aorta <br>
Gene1	SNP1	-0.0350153	-0.0796675	0.0458593	-0.0663155 <br>
Gene2	SNP2	0.25088	0.133673	0.13425	0.211878 <br>
Gene3	SNP3	0.0262571	-0.065221	0.199401	-0.0711795 <br>
Gene4	SNP4	-0.272452	0.240933	0.214758	0.281942 <br>
Gene5	SNP5	NA	NA	NA	NA <br>
Gene6	SNP6	0.133723	0.0933188	0.103415	-0.15649 <br>
```


```data/test_data_W.txt```: each row contains the weight (reciprical of standard error of the effect size) of an eQTL across tissues. Columns should be aligned with the columns in ```data/test_data_X.txt```.

```
Gene	SNP	Adipose_Subcutaneous	Adipose_Visceral_Omentum	Adrenal_Gland	Artery_Aorta <br>
Gene1 SNP1	0.0748711	0.0926145	0.150558	0.0754927 <br>
Gene1	SNP2	0.0425708	0.036122	0.0405176	0.0548538 <br>
Gene1	SNP3	0.0735933	0.0765909	0.125968	0.0891406 <br>
Gene1	SNP4	0.164811	0.152243	0.235161	0.177724 <br>
Gene1	SNP5	NA	NA	NA	NA <br>
Gene1	SNP6	0.114314	0.112615	0.182777	0.147263 <br>
```


### Model selection


Model selection is one of the most challenging parts in deciding matrix factorization models. People have used several methods to approach this problem (REF: xxxxxx). In sn-spMF, we recommend searching for the hyper-parameters (K, alpha1, lambda1) in two steps:

Narrow down the range of hyper-parameters 

When first running the algorithm, it may be completely unclear how to choose the appropriate range to search for hyper-parameters. We recommend first searching for the appropriate range, by 1). running the scripts in well-separated numerical ranges, like choose from [1, 10, 100, 500, 1000]; and 2).  setting the number of iterations to a moderate number since there is no need to reach accurate results. 

If the number of factors become much smaller than the initial number of factors to start with (ie. a lot of factors become empty), it means that the penalty parameters are too stringent. Usually we have an estimated level of sparsity, for example, around 80%, for the loading matrix and factor matrix. If the reported sparsity is far below the expected sparsity (ie. 20%), it means that the penalty parameters are too small. 

Based on the initial round of searching, we should have the numerical range to search for. 


Refine the hyper-parameter selection

With the learned range of hyper-parameters, we continue to look in finer grids. For example, run the scripts for alpha1 in [10, 20, 30, … 100]. Because the three parameters can collaboratively affect the decomposition results, we perform model selection in two sub-steps:

2.1). Choose the number of factors K. 

We notice that the cophenetic coefficient can be affected by sparsity in the decomposed matrices given different settings of alpha1 and lambda1 with fixed K. To gain more stable matrix decomposition results, we compare the average cophenetic coefficient with multiple settings for alpha1 and lambda1. The optimal K is chosen to have the highest average cophenetic coefficient. 

2.2). Choose the penalty parameters alpha1 and lambda 1. 

Because factors are expected to be independent of each other, to alleviate multicollinearity, we then search for the alpha1 and lambda1 that result in factors with smallest correlation. 


We’d like to include some suggestions from practical experience when setting the arguments in the model:

Number of iterations: recommend using 100 or less. If the model does not converge within 100 iterations, one reason is that the penalty parameters are too small, which leads to slow optimization steps. Larger penalty parameters are suggested in the case where the model does not converge within 100 iterations. 

Change in the factor matrix to call convergence (converged_F_change): this is the Frobenius norm of the difference matrix comparing the factor matrix before and after updating, scaled by the number of factors (||F_new - F_old||^2_F / (number of factors)). The scaling is to avoid bias of higher Frobenius norm coming from more factors. 

Change in the objective to call convergence (converged_obj_change): this is usually a more stringent threshold than converged_F_change. 

Number of runs to compute cophenetic coefficient: we find that around 20 runs suffice to provide a reliable estimate of the cophenetic coefficient. 


