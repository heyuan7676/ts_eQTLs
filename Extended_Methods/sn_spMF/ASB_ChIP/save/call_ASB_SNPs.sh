#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


ml python/2.7
python call_ASB_SNPs.py "$1" "$2"
