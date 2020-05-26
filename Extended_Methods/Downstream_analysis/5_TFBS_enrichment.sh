#!/bin/bash
#SBATCH --time 40:00:00
#SBATCH --no-requeue
#SBATCH -p skylake
#SBATCH --nodes=1
#SBATCH --ntasks=1

source ./GLOBAL_VAR.sh

module load python/2.7
python ${scripts_dir}/5_TFBS_enrichment.py "$1" "$2" "$3"
