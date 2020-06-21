#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

ml python/2.7
python 5_TFBS_ChIP_ASB-HM.py "$1"









