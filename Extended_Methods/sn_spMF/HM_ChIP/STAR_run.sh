#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

bash STAR_bamfile_format.sh
bash STAR_intersect_outlier_SNPs.sh
