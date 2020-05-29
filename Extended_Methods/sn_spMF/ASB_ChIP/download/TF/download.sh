#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

xargs -L 1 curl -O -L < download_filtered.txt
