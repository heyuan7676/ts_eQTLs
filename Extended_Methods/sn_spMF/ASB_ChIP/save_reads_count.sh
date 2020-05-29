#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p shared


ml python/2.7
python save_reads_count.py "$1" "$2" 1
