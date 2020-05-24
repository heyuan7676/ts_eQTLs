#!/bin/bash
#SBATCH --time 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p skylake


### Search the appropriate range for alpha1 and lambda1
### 1). search in well-seperated ranges
### 2). set iterations to a moderate number, there is no need to reach the accurate results

K="$1"
alpha1="$2"
lambda1="$3"

iterations="$4"


ml R/3.5.1
Rscript sn_spMF/run_MF.R -k $K -a ${alpha1} -l ${lambda1} -t ${iterations}


### Examine the printed output to narrow down the proper range of alpha1 and lambda1 for further search
