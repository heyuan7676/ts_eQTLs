This folder contains scripts used to generate simulations and to fit the simulation with different models. 
Breif description of the scripts are as below. Please find the details inside the scripts.

##### Note: all scripts should be run from the root dir. (ie. ts_eQTLs/)

### Prerequisites
R code is run in ```R/3.5.1```. 

R packages needed are:
```
install.packages('NMF')
install.packages('PMA')
install.packages('ashr')
install.packages('flashr')
install.packages('softImpute')
install.packages('ssvd')
install.packages('combinat')
install.packages('R.matlab')
install.packages('RColorBrewer')
install.packages('cowplot')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('pheatmap')
install.packages('plyr')
install.packages('reshape2')
install.packages('optparse')
```


```Generate_input.R```: Generate input data in the simulation. 

```choose_hyperparameters_simulate.sh```: Run sn_spMF to fit the simulated data, perform model selection.

```tune_parameters.R```: Wrapper to tune parameters. Calls ```sn_spMF/collect_results.R```

```tune_parameters_function.R```: Wrapper to choose the optimal hyper-parameters and optimal solution. 

```perform_MF_methods.R```: Run various models. For sn_spMF, ```tune_parameters_function.R``` is called. 

```fit_significant_hits.R```: Map eQTLs to the learned factors.

```compare_methods.R```: Wrapper to compare results from various models. 

```compare_methods_simulate.sh```: Call ```compare_methods.R```.


### Run simulation and evaluate 

#### 1. Generate the input files
```
## generate the input
Rscript simulation/Generate_input.R -t ${tau} -s ${seed}
```

#### 2. Tune parameters for sn_spMF
```
## tune parameters for sn_spMF
sbatch simulation/choose_hyperparameters_simulate.sh ${tau} ${seed}

## After all choose_hyperparameters_simulate.sh finish running -
Rscript simulation/tune_parameters.R -O simulation/output/tau${tau}_seed${seed}/
```

#### 3. Run all models, and collect the metrics
```
## collect results from all methods
bash simulation/compare_methods_simulate.sh ${tau} ${seed}
```

#### 4. Visualize the metrics
```
## plot the resulting metrics
Rscript simulation/plot_metrics.R 

```
