#!/bin/bash
#SBATCH --time 30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue

source ./GLOBAL_VAR.sh
NBlocks=38
ml python/2.7
ml parallel


### Step 1)
### aggregate pairs for credible sets


### Step 2.2)
### collect cis-eQTL results for all eQTL pairs in the credible set
### parallel implementation in the script
bash step2_collect_cisresults_bg.sh


### Step 3)
### Merge info across tissues
bash step3_merge_files_bg.sh


