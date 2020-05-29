#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --no-requeue


xargs -n 1 curl -O -L < download_filtered.txt



## quickly check the metadata.tsv
#for f in *bed.gz; do x=${f%.bed.gz}; cat metadata.tsv | grep "^${x}" | wc -l; done

