This folder contains scripts used to generate simulations and to fit the simulation with different models. 
Breif description of the scripts are as below. Please find the details inside the scripts. <\br>

##### Note: all scripts should be run from the root dir. (ie. ts_eQTLs/)


```Generate_input.R```: Generate input data in the simulation. 

```choose_hyperparameters_simulate.sh```: Run sn_spMF to fit the simulated data, perform model selection.

```tune_parameters.R```: Wrapper to tune parameters. Calls ```sn_spMF/collect_results.R```

```tune_parameters_function.R```: Wrapper to choose the optimal hyper-parameters and optimal solution. 

```perform_MF_methods.R```: Run various models. For sn_spMF, ```tune_parameters_function.R``` is called. 

```fit_significant_hits.R```: Map eQTLs to the learned factors.

```compare_methods.R```: Wrapper to compare results from various models. 

```compare_methods_simulate.sh```: Call ```compare_methods.R```.


To run simulation for a specific ```tau``` and ```seed```, please run:
```
## generate the input
Rscript simulation/Generate_input.R -t ${tau} -s ${seed}

## tune parameters for sn_spMF
sbatch simulation/choose_hyperparameters_simulate.sh ${tau} ${seed}

## After all choose_hyperparameters_simulate.sh finish running -
Rscript simulation/tune_parameters.R -O simulation/output/tau${tau}_seed${seed}/

## collect results from all methods
bash simulation/compare_methods_simulate.sh ${tau} ${seed}

```
