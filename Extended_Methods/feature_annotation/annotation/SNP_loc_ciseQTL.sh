#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=2:00:00
#SBATCH -p debug

source ./GLOBAL_VAR.sh
datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets

inputdir=${datadir}/cbset_datasets/caviar_output_GTEx_LD/aggregate
outputdir=${inputdir}/SNP_loc



cd ${inputdir}
rm -f -r ${outputdir}
mkdir -p ${outputdir}
for f in `ls *_95set_pairs.txt`
do
	echo $f
	outfn=SNP_loc_${f%.txt}.bed
	paste <(cat ${f} | awk '{print $2}' | awk -F'_' '{print $1}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=2 '{print $2-s}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=2 '{print $2+s}')  <(cat ${f} | awk '{print $2}') >> ${outputdir}/${outfn}
done



cd ${outputdir}
fn=commands.txt
for f in `ls SNP_loc_*bed`
do
	echo "sort -k1,1 -k2,2n ${f} | uniq > temp_${f}" >> ${fn}
done

ml parallel
parallel :::: ${fn}


for f in `ls temp_*`
do
	mv ${f} ${f#temp_}
done	

