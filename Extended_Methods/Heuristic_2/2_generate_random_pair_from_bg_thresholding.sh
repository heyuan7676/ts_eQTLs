#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load python/2.7
python 2_generate_random_pair_from_bg_thresholding.py "$1" "$2"
