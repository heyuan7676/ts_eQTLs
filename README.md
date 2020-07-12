
# Learn latent factors using sn-spMF
We develop a constrained matrix factorization model to learn patterns of tissue-sharing and tissue-specificity of eQTLs across 49 human tissues from the Genotype-Tissue Expression (GTEx) project. The learned factors include patterns reflecting tissues with known biological similarity or shared cell types, in addition to a dense factor representing a ubiquitous genetic effect across all tissues.

## Prerequisites
R code is run in ```R/3.5.1```. 

Install required R packages:
```
pkgs.list <- c("penalized", "readr", "plyr", "reshape2", "optparse", 
                "dplyr", "colorspace", "colortools", "cowplot", "ggplot2",
                "gridExtra", "lemon")
new.pkgs <- pkgs.list[!(pkgs.list %in% installed.packages()[,"Package"])]
if(length(new.pkgs) > 0 ) { 
	install.packages(new.pkgs)
}
```

## Run the sn_spMF model
To get the result for one run, please run the following command. Details can be found in ```sn_spMF/run_MF.R```.
```
Rscript sn_spMF/run_MF.R -k 17 -a 100 -l 50 -t 100
```

### Input
There are two important features of input files for sn_spMF:

##### 1). sn_spMF is able to learn the underlying patterns from subset of data
For example, lead eQTLs among all eQTLs in the credible set. To demonstrate this, we provide the demo data as in ```data/test_data_X_all.txt``` and ```data/test_data_SE_all.txt```, and derived a subset of all the eQTLs as in ```data/test_data_X.txt``` and ```data/test_data_SE.txt```. We used the subset of data points to learn the factor matrix, and then map all eQTLs to the factors.

##### 2). Allow missing data in the input. 
For each data point (or eQTL), it can have missing data in various tissues. Removing data points with any missing data will not only loss information, but also cause bias by focusing more on the shared eQTLs. sn_spMF takes care of missing data by assigning weights of zero when computing the objective, and thus avoid removing missing data. 


##### Demo of the input files.

```data/test_data_X.txt```: each row contains the effect size of an eQTL across tissues; the first two columns are gene names and SNP names for the eQTLs, and following columns are the features to learn patterns about, (tissues in the demo, can be time points in time-series data, or cells in single cell data). Missing data are presented as NA. Columns are separated by '\t'. 

```
Gene	SNP	Adipose_Subcutaneous	Adipose_Visceral_Omentum	Adrenal_Gland	Artery_Aorta
Gene1	SNP1	-0.0350153	-0.0796675	0.0458593	-0.0663155
Gene2	SNP2	0.25088	0.133673	0.13425	0.211878
Gene3	SNP3	0.0262571	-0.065221	0.199401	-0.0711795
Gene4	SNP4	-0.272452	0.240933	0.214758	0.281942
Gene5	SNP5	NA	NA	NA	NA
Gene6	SNP6	0.133723	0.0933188	0.103415	-0.15649
```


```data/test_data_SE.txt```: each row contains the standard error of the effect size of an eQTL across tissues. Columns should be aligned with the columns in ```data/test_data_X.txt```.

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

Users can find the learned factor matrix in output/sn_spMF_K17_a1100_l150/sn_spMF_K17_a1100_l150.\* including the plotted factors. The output dir can be specified using ```-O``` when running  ```sn_spMF/run_MF.R```.


## (Optional) Multiple initializations

Because random initializations can result in different decomposition solutions, we recommend running the decomposition multiple times (ie. 30 times), and obtain the optimal solution using the decomposition with minimum objective value. Users can run the following to extract the solution with optimal objective (saved in ```output/sn_spMF_K17_a1100_l150/*RData``` by default, can be changed using the ```-O``` argument), or extract the solution with optimal objective from the model selection step (see below, section "Model Selection").

```
## Run intialization multiple times
Rscript sn_spMF/run_MF.R -k 17 -a 100 -l 50 -t 100 -c 1

## Extract the optimal solution
Rscript sn_spMF/find_optimal.R -k 17 -a 100 -l 50 
```

The resulting optimal solution for factor matrix in this implementation looks like:

![alt text](https://github.com/heyuan7676/ts_eQTLs/blob/master/output/sn_spMF_K17_a1100_l150/sn_spMF_FactorMatrix_K17_a1100_l150_Run25_with_legends.png)
## Model selection

In the sn-spMF model, we need to set hyper-parameters including the rank of the decomposition (K) and the sparsity penalties (alpha, lambda). We recommend searching for the hyper-parameters (K, alpha, lambda) in two steps:

#### 1. Narrow the sparsity penalty hyper-parameter search space

Exploring three hyper-parameters jointly can be computationally expensive, and many values would produce clearly implausible models. In order to choose a tractable, appropriate range of hyper-parameter settings to evaluate, we recommend running the method for a broad range of well-separated settings of sparsity penalty hyper-parameters, such as  [1, 10, 100, 500], and a wide range of K chosen by considering the number of tissues or experiments in the data. This coarse-grained search step can be evaluated using the following guidelines:

a). the sparsity of the solutions being in accordance with user expectations for their domain - if the reported sparsity is far below the expected sparsity, the chosen penalty parameters may be too small.

b). behavior of the factor matrix: if the number of utilized factors become much smaller than the initial number of factors to start with (ie. a lot of factors become empty, having no non-zero entries), it means that the penalty parameters are likely too stringent. 
An example to perform this step is as below:
```
iterations=20
for K in 10 15 20
do
        for alpha in 1 10 100 500
        do
                for lambda in 1 10 100 500
                do
                        bash sn_spMF/1_run_parameter_scope_search.sh ${K} ${alpha} ${lambda} ${iterations}
                done
        done
done
```

To collect the results from multiple runs, users can run the following command. The output will be saved in output/choose_para_preliminary.txt
```Rscript sn_spMF/tune_parameters_preliminary.R -f choose_para_preliminary.txt```

We ran the code above for demo data. When examining the output file output/choose_para_preliminary.txt, we observe that  and  of either 1 or 10 result in a factor matrix with sparsity of 10% - 50%, which is lower than our expectation for this example (or a multi-tissue domain such as GTEx). On the other hand, alpha and lambda of 500 result in too few utilized, non-zero factors (for example, K=10, alpha=500, lambda=1 result in only around 6 non-zero factors).   and  of 10 or 100 appear to give a balance between sparsity in the factor matrix and number of non-zero factors. K=10 results in the majority of solutions having 10 used factors, while K=20 results in the majority of solutions having nearly 20 used factors, after eliminating the models with sparsity penalties that are too stringent. This shows that searching between 10 and 20 for K is reasonable in this situation. Thus we proceed to perform grid search for alpha and lambda in the range of 10 to 100, and K in the range of 10 to 20.


#### 2. Refine the sparsity penalty hyper-parameter selection

Within a manageable search space for the hyper-parameters as selected above, we then suggest searching settings using finer granularity and evaluating the learned models for stability along with independence between factors. For example, run for alpha and lambda in [10, 20, 30, ... 100] or finer, and run the model multiple times from random initializations (ie. 30 times). We recommend using a maximum number of iterations for 100 or less for each run. If the model does not converge within 100 iterations, it is probably because the penalty parameters are too small, which leads to very slow optimization steps. Larger penalty parameters are suggested in the case where the model does not converge within 100 iterations. An example to perform this step is as below. 
```
iterations=100
for K in {10..20}
do
        for alpha in {1..10}
        do
                for lambda in {1..10}
                do
                        a=$(( 10*alpha ))
                        l=$(( 10*lambda ))
                        run sn_spMF/2_choose_hyperparameters.sh ${K} ${a} ${l} ${iterations}
                done
        done
done
```

Within these chosen search spaces, we evaluated sn-spMF models for all combinations of K, alpha and lambda using 1) a previously defined criterion of matrix factorization stability by Brunet et al. [1], and 2) independence of the learned factors, which represents adequate sparsity. Considering the stochastic nature of matrix factorization, Brunet et al. proposed a method looking for the most stable factorization result, and this method has been applied in various studies [1,2]. We obtained the consensus matrix C after 30 runs with random initialization for each model. The values in C are between 0 to 1, representing the proportion of runs in which a pair of tissues are assigned to the same factor. Using the C matrix, we computed the cophenetic correlation which is used to measure the degree of dispersion for the C matrix. Higher cophenetic correlation indicates a more stable factor matrix. To collect the evaluation metrics, users can run the following command. The output will be saved in output/choose_para.txt.

```Rscript sn_spMF/tune_parameters.R -f choose_para.txt```

Based on the evaluation metrics, we performed the following selection steps: 

a). We first eliminated some settings of K.  Here, for each observed mean number of learned, non-empty factors K' (which may be less than the input K), we aggregated across the different settings of alpha and lambda and computed the median cophenetic correlation [1].  

b). We eliminated from consideration any settings of K corresponding to a K' with a median cophenetic correlation <0.9. Next, among the remaining individual settings, we eliminated any cophenetic correlation <0.9.  

c). Last, among these apparently stable settings, we selected the final hyper-parameters based on the minimum Pearson correlation between pairs of factors, to encourage independent factors and a level of sparsity that matches independent signals in the data. Here, we computed Pearson correlation for each pair of factors, took the Frobenius norm of the pairwise correlation matrix, and averaged this across the 30 randomly initialized runs for the same setting.   

In the demo data, We chose settings of K corresponding to a K' higher than 9, such that the corresponding median cophenetic correlation is above 0.9, and followed steps b) and c) to select the optimal model solution. The script is available in sn_spMF/choose_paras_sn_spMF.ipynb. A separate example of learning the hyper-parameters is provided in simulation/choose_paras_sn_spMF_simulation.ipynb on simulated data. Details can be found in simulation/.


![alt text](https://github.com/heyuan7676/ts_eQTLs/blob/master/output/choose_para_K.png)
![alt text](https://github.com/heyuan7676/ts_eQTLs/blob/master/output/choose_para_nFactors.png)


## Examine the optimal solution.

By examining the tuning results in ```sn_spMF/choose_paras_sn_spMF.ipynb```, we find that ```sn_spMF_FactorMatrix_K17_a1100_l150``` is the optimal setting of hyper-parameters. Among the 30 runs using this implementation, ```run25``` gives the optimal solution with the minimum objective. Users can find the learned factor matrix in ``` output/sn_spMF_K17_a1100_l150/sn_spMF_K17_a1100_l150_Run25.*```, including the plotted factors. 

The resulting optimal solution for factor matrix looks like:

![alt text](https://github.com/heyuan7676/ts_eQTLs/blob/master/output/sn_spMF_K17_a1100_l150/sn_spMF_FactorMatrix_K17_a1100_l150_Run25_factors.png)




## Map eQTLs to factors.
After user have chosen the optimal hyper-parameters (```${FM_fn}```), please run the following command to map the eQTLs to the learned factors. The script automatically chooses the solution with the optimal objective if multiple solutions exist. The mapped eQTLs are in ```output/mapping/``` by default or can be specified by ```-m ${mappingDir}```. Details can be found in ```mapping/lm.R```.

```
K=17
alpha1=100
lambda1=50
FM_fn=sn_spMF_K${K}_a1${alpha1}_l1${lambda1}
Rscript mapping/lm.R -f ${FM_fn}
```

## Brief description of folders in this repository
`sn_spMF/`: main folder with code for running sn-spMF to learn latent patterns.

`mapping/`: map eQTLs to factors after learning the latent patterns. 

`data/`: demo data

`output/`: output from running experiment on demo data, including inter-mediate results for paramter selection.

`simulation/`: run different matrix factorization methods on simulated data.

`Extended_Methods/`: code used in the paper, including heuristic methods, and downtream analysis

`plots_in_the_paper/`: code used to generate figures in the paper


## Reference
[1]. Brunet, J.-P., Tamayo, P., R Golub, T., P Mesirov, J.: Metagenes and
molecular pattern discovery using matrix factorization. Proceedings of
the National Academy of Sciences. 101, 4164–9 (2004). doi:10.1073/pnas.0308531101

[2]. Wu, S., Joseph, A., S. Hammonds, A., E. Celniker, S., Yu, B., Frise,
E.: Stability-driven nonnegative matrix factorization to interpret spatial
gene expression and build local gene networks. Proceedings of the
National Academy of Sciences 113, 201521171 (2016). doi:10.1073/pnas.1521171113


