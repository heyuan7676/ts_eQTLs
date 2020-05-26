This folder contains scripts to run sn_spMF. 
Breif description of the scripts are as below. Please find the details inside the scripts.

##### Note: all scripts should be run from the root dir. (ie. ts_eQTLs/)

```run_MF.R```: Wrapper to run sn_spMF. User can specify the hyper-parameters, input files, output dir, as well as convergence criteria to run sn_spMF. Some arguments include:

        - ```Number of iterations```: we recommend using 100 or less. If the model does not converge within 100 iterations, one reason is that the penalty parameters are too small, which leads to slow optimization steps. Larger penalty parameters are suggested in the case where the model does not converge within 100 iterations. 

        - ```Change in the factor matrix to call convergence (converged_F_change)```: this is the Frobenius norm of the difference matrix comparing the factor matrix before and after updating, scaled by the number of factors (||F_new - F_old||^2_F / (number of factors)). The scaling is to avoid bias of higher Frobenius norm coming from implementations with more factors. 

        - ```Change in the objective to call convergence (converged_obj_change)```: this is usually a more stringent threshold than converged_F_change. 

        - ```Number of runs to compute cophenetic coefficient```: we find that around 20-30 runs suffice to provide a reliable estimate of the cophenetic coefficient. 


```sn_spMF.R```: Main function to run weighted sn-spMF. 


```readIn.R```: Read in the input files.

```Update_FL.R```: Wrapper to run optimization in sn_spMF. Main functions include ```fit_L.R``` and ```fit_F.R```. 

```cophenet.R``` and ```compute_obj.R```: Helper functions to compute the objective, and the cophenetic coefficient. 

```plot_input_matrix.R``` and ```plot_factor_matrix.R```: Helper functions to visualize the input matrix and learned factor matrix.

```collect_results.R```: Function to collect the evaluation metrics for parameter tuning. 

```tune_parameters_preliminary.R``` and ```tune_parameters.R```: Wrapper to tune parameters. Calls ```collect_results.R```

```find_optimal.R```: Given a setting of hyper-parameters, extract the optimal run with minimum objective
