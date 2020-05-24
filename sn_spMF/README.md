This folder contains scripts to run sn_spMF.<br>

```run_MF.R```: Wrapper to run sn_spMF. User can specify the hyper-parameters, input files, output dir, as well as convergence criteria to run sn_spMF. <br>
Example of usage:
```
Rscript run_MF.R -k 20 -a 10 -l 200 -t 10 -c 1 -r 1
```

```sn_spMF.R```: Main function to run weighted sn-spMF. Details about input can be found in the script. 


```readIn.R```: Read in the input files.

```Update_FL.R```: Wrapper to run optimization in sn_spMF. Main functions include ```fit_L.R``` and ```fit_F.R```. 

```cophenet.R``` and ```compute_obj.R```: Helper functions to compute the objective, and the cophenetic coefficient. 

```plot_input_matrix.R``` and ```plot_factor_matrix.R```: Helper functions to visualize the input matrix and learned factor matrix.

```collect_results.R```: Function to collect the evaluation metrics for parameter tuning. 

```tune_parameters_preliminary.R``` and ```tune_parameters.R```: Wrapper to tune parameters. Calls ```collect_results.R```
