#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=2

module load python/2.7
python 5_TFBS_snps_genes_thr.py "$1" "$2"
