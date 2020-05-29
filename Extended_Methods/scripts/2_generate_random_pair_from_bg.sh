#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1


source ./GLOBAL_VAR.sh

module load python/2.7
python ${scripts_dir}/2_generate_random_pair_from_bg.py "$1" "$2"
