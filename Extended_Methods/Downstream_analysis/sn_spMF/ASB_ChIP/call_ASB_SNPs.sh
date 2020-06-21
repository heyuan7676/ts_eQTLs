#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p shared

ml python/2.7
python call_ASB_SNPs.py HNF4A 10 2
python call_ASB_SNPs.py CTCF 10 2
