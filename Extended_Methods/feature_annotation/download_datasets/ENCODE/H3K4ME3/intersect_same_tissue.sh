#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

module load bedtools

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/ENCODE
cd ${datadir}
rm -f prefix.txt

for x in *bed
do
        echo ${x} | cut -d'_' -f1,2 >> prefix.txt
done

cat prefix.txt | sort | uniq > temp
mv temp prefix.txt


while IFS= read -r line
do
        t=$line
        echo $t
	rm -f ${t}_union.bed
	rm -f ${t}_intersect.bed
        cat ${t}*bed | sort -k1,1 -k2,2n | bedtools merge -c 5 -o median > ${t}_union.bed
	multiIntersectBed -i ${t}*bed > ${t}_intersect.bed
done < prefix.txt
