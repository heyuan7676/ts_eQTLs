#!/bin/bash
#SBATCH --time 40:00:00
#SBATCH --no-requeue
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=2

source ./GLOBAL_VAR.sh

module load python/2.7
python ${scripts_dir}/5_TFBS_enrichment_results.py "$1" "$2" "$3"
