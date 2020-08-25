#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

ml python/2.7
python save_ASB.py "$1" "$2" 2
