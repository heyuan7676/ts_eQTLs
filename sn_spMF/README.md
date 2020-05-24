This folder contains scripts to run sn_spMF. Breif description of the scripts are as below. Please find the details inside the scripts.

```run_MF.R```: Wrapper to run sn_spMF. User can specify the hyper-parameters, input files, output dir, as well as convergence criteria to run sn_spMF.

```sn_spMF.R```: Main function to run weighted sn-spMF. 


```readIn.R```: Read in the input files.

```Update_FL.R```: Wrapper to run optimization in sn_spMF. Main functions include ```fit_L.R``` and ```fit_F.R```. 

```cophenet.R``` and ```compute_obj.R```: Helper functions to compute the objective, and the cophenetic coefficient. 

```plot_input_matrix.R``` and ```plot_factor_matrix.R```: Helper functions to visualize the input matrix and learned factor matrix.

```collect_results.R```: Function to collect the evaluation metrics for parameter tuning. 

```tune_parameters_preliminary.R``` and ```tune_parameters.R```: Wrapper to tune parameters. Calls ```collect_results.R```
