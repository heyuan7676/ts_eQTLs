
# Learn latent factors using sn-spMF
We develop a constrained matrix factorization model to learn patterns of tissue-sharing and tissue-specificity of eQTLs across 49 human tissues from the Genotype-Tissue Expression (GTEx) project. The learned factors include patterns reflecting tissues with known biological similarity or shared cell types, in addition to a dense factor representing a ubiquitous genetic effect across all tissues.

## Prerequisites
R code is run in ```R/3.5.1```. 

R packages needed are:
```
install.packages('penalized')
install.packages('readr')
install.packages('plyr')
install.packages('reshape2')
install.packages('optparse')
install.packages('dplyr')
install.packages('colorspace')
install.packages('colortools')
install.packages('cowplot')
install.packages('ggplot2')
install.packages('gridExtra')
install.packages('lemon')
```

## Run the sn_spMF model
To get the result for one run, please run the following command. Details can be found in ```run_MF.R```.
```
Rscript sn_spMF/run_MF.R -k 17 -a 100 -l 50 -t 100
```

### Input

sn_spMF is able to learn the underlying patterns from subset of data, for example, lead eQTLs among all eQTLs in the credible set. To demonstrate this, we provide the demo data as in ```data/test_data_X_all.txt``` and ```data/test_data_SE_all.txt```, and derived a subset of all the eQTLs as in ```data/test_data_X.txt``` and ```data/test_data_SE.txt```. We used the subset of datapoints to learn the factor matrix, and then map all eQTLs to the factors.

```data/test_data_X.txt```: each row contains the effect size of an eQTL across tissues; the first two columns are gene names and SNP names for the eQTLs, and following columns are the features to learn patterns about, (tissues in the demo, can be time points in time-series data, or cells in single cell data). Missing data are presented as NA. Columns are seperated by '\t'. Model sn_spMF takes care of missing data by assigning weights of zero when computing the objective.

```
Gene	SNP	Adipose_Subcutaneous	Adipose_Visceral_Omentum	Adrenal_Gland	Artery_Aorta
Gene1	SNP1	-0.0350153	-0.0796675	0.0458593	-0.0663155
Gene2	SNP2	0.25088	0.133673	0.13425	0.211878
Gene3	SNP3	0.0262571	-0.065221	0.199401	-0.0711795
Gene4	SNP4	-0.272452	0.240933	0.214758	0.281942
Gene5	SNP5	NA	NA	NA	NA
Gene6	SNP6	0.133723	0.0933188	0.103415	-0.15649
```


```data/test_data_W.txt```: each row contains the weight (reciprical of standard error of the effect size) of an eQTL across tissues. Columns should be aligned with the columns in ```data/test_data_X.txt```.

```
Gene	SNP	Adipose_Subcutaneous	Adipose_Visceral_Omentum	Adrenal_Gland	Artery_Aorta
Gene1 SNP1	0.0748711	0.0926145	0.150558	0.0754927
Gene1	SNP2	0.0425708	0.036122	0.0405176	0.0548538
Gene1	SNP3	0.0735933	0.0765909	0.125968	0.0891406
Gene1	SNP4	0.164811	0.152243	0.235161	0.177724
Gene1	SNP5	NA	NA	NA	NA
Gene1	SNP6	0.114314	0.112615	0.182777	0.147263
```

### Output

User can find the learned factor matrix in output/sn_spMF_K17_a1100_l150/sn_spMF_K17_a1100_l150.\* including the plotted factors. The output dir can be specificed using ```-O``` when running  ```sn_spMF/run_MF.R```.


## (Optional) Multiple intializations

Because random initializations can result in different decomposition solutions, we recommend running the decomposition multiple times (ie. 30 times), and obtain the optimal solution using the decomposition with minimum objective value. User can run the following to extract the solution with optimal objective (saved in ```output/sn_spMF_K17_a1100_l150/*RData``` by default, can be changed using the ```-O``` argument), or extract the solution with optimal objective from the model selection step (see below).

```
## Run intialization multiple times
Rscript sn_spMF/run_MF.R -k 17 -a 100 -l 50 -t 100 -c 1

## Extract the optimal solution
Rscript sn_spMF/find_optimal.R -k 17 -a 100 -l 50 
```

The resulting factor solution looks like:

![alt text](https://github.com/heyuan7676/ts_eQTLs/blob/master/output/sn_spMF_K17_a1100_l150/sn_spMF_FactorMatrix_K17_a1100_l150_Run25_factors.png)


## Model selection


Model selection is one of challenging steps in constructing matrix factorization models. In sn-spMF, we recommend searching for the hyper-parameters (K, alpha1, lambda1) in two steps:

#### 1. Narrow down the range of hyper-parameters 

When first running the algorithm, it may be completely unclear how to choose the appropriate range to search for hyper-parameters. We recommend first searching for the appropriate range, by 1). running the scripts in well-separated numerical ranges, like choose from [1, 10, 100, 500, 1000]; and 2).  setting the number of iterations to a moderate number since there is no need to reach accurate results. 

If the number of factors become much smaller than the initial number of factors to start with (ie. a lot of factors become empty), it means that the penalty parameters are too stringent. Usually we have an estimated level of sparsity, for example, around 80%, for the loading matrix and factor matrix. If the reported sparsity is far below the expected sparsity (ie. 20%), it means that the penalty parameters are too small. 

Based on the initial round of searching, we should have the numerical range to search for. 

An example to perform this step is as below:
```
iterations=20
for K in 10 15 20
do
        for alpha1 in 1 10 100 500
        do
                for lambda1 in 1 10 100 500
                do
                        sbatch 1_run_parameter_scope_search.sh ${K} ${alpha1} ${lambda1} ${iterations}
                done
        done
done
```

To collect the results from multiple runs, user can run the following command. The output will be saved in ```output/choose_para_preliminary.txt```
```
Rscript sn_spMF/tune_parameters_preliminary.R -f choose_para_preliminary.txt
```



#### 2. Refine the hyper-parameter selection

With the learned range of hyper-parameters, we continue to look in finer grids. For example, run the scripts for alpha1 and lambda1 in [10, 20, 30, … 100]. An example to perform this step is as below. Detailed description of arguments in the model are described in sn_spMF/.
```
iterations=100
for K in {10..20}
do
        for alpha1 in {1..10}
        do
                for lambda1 in {1..10}
                do
                        a=$(( 10*alpha1 ))
                        l=$(( 10*lambda1 ))
                        sbatch 2_choose_hyperparameters.sh ${K} ${a} ${l} ${iterations}
                done
        done
done
```

To collect the results from multiple runs, user can run the following command. The output will be saved in ```output/choose_para.txt```.
```
Rscript sn_spMF/tune_parameters.R -f choose_para.txt
```

Because the three parameters can collaboratively affect the decomposition results, we perform model selection in two sub-steps, for which we provide an example in ```choose_paras_sn_spMF.ipynb```. 


##### 2.1 Choose the range of number of factors. 

We notice that the cophenetic coefficient can be affected by sparsity in the decomposed matrices given different settings of alpha1 and lambda1 with fixed K. To gain more stable matrix decomposition results, we compare the average cophenetic coefficient with multiple settings for alpha1 and lambda1. In the demo data, different implementations of K doesn't result in obvious difference in the cophenetic coefficient, and thus we do not to filter on K. We observe that some implementations push factors to be zero and thus the real number of factors reached is different from the assigned number of factors. We choose number of learned factors to be those with median cophenetic coefficient > 0.9. 

##### 2.2 Filter out implementations with low cophenetic coefficient

We then filter out implementations with cophenetic coefficient < 0.9, and keep implementations with consistent decomposition solutions given random initializations indicating the stability of the solution. In this situation, we try different values of the thresholds including 0.85 and 0.8, and observe that the same optimal solution is learned. 

##### 2.3 Choose the implementation with the most independent factors

Because factors are expected to be independent of each other, to alleviate multicollinearity, we then search for the alpha1 and lambda1 that result in factors with smallest correlation. 

##### Note: a seperate example of learning the hyper-parameters is provided in ```simulation/choose_paras_sn_spMF_simulation.ipynb``` on simulated data. Details can be found in ```simulation/```. 

## Examine the optimal solution.

By examining the tuning results in ```choose_paras_sn_spMF.ipynb```, we find that ```sn_spMF_FactorMatrix_K17_a1100_l150``` is the optimal setting of hyper-parameters. Among the 30 runs using this implementation, ```run25``` gives the optimal solution with the minimum objective. User can find the learned factor matrix in ``` output/sn_spMF_K17_a1100_l150/sn_spMF_K17_a1100_l150_Run25.*```, including the plotted factors. 



## Map eQTLs to factors.
After user have chosen the optimal hyper-parameters (```${FM_fn}```), please run the following command to map the eQTLs to the learned factors. The script automatically chose the solution with optimal objective if multiple solutions exist. The mapped eQTLs are in ```output/mapping/``` by default or can be specified by ```-m ${mappingDir}```. Details can be found in ```mapping/lm.R```.

```
K=17
alpha1=100
lambda1=50
FM_fn=sn_spMF_K${K}_a1${alpha1}_l1${lambda1}
Rscript mapping/lm.R -f ${FM_fn}
```

