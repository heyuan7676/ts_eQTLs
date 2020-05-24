This folder contains scripts used to map eQTLs to the learned factors.
Breif description of the scripts are as below. Please find the details inside the scripts.

Note: all scripts should be run from the root dir. (ie. ts_eQTLs/)


```readin_data.R```: Wrapper to read in the input files, factor matrix and learned loading matrix.

```run_linearReg.R```: Perform weighted linear regression. 

```adj_pvalue.R```: Perform Benjamini-Hochberg Procedure to control for the false discover rate. 

```Remove_compensate_pairs.R```: Remove eQTLs with opposite loadings compared to the effect size. Details can be found in section Methods - "Assignment of eQTLs to factors"
in the paper. 


```lm.R```: Wrapper to map eQTLs to the factors. 
