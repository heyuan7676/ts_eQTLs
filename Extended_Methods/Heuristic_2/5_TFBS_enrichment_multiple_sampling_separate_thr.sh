#!/bin/bash
#SBATCH --time 40:00:00
#SBATCH --no-requeue
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load python/2.7

python 5_TFBS_enrichment_multiple_sampling_separate_noRes_thr.py "$1" "$2" "$3"
