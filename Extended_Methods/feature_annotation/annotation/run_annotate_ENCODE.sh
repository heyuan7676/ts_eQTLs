#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00

hm="$1"
for chr in {1..22}
do 
	echo $chr
	bash annotate_ENCODE.sh $hm union $chr
	bash annotate_ENCODE.sh $hm intersect $chr
done
