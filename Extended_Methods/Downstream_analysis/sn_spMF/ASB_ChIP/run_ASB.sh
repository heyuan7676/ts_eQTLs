#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p shared

ml python/2.7

TF="$1"
thr="$2"
peaks="$3"

python save_ASB.py ${TF} ${thr} ${peaks}
python call_ASB_SNPs.py ${TF} ${thr} ${peaks}
