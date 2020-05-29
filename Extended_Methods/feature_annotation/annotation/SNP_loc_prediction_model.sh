#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=2:00:00
#SBATCH -p debug

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets

inputdir=${datadir}/GTEx_datasets/allPairs/maf_dis/split_by_gene/allpairs_including_outliers
outputdir=${datadir}/prediction_model/SNP_loc



cd ${inputdir}
rm -r -f ${outputdir}
mkdir -p ${outputdir}
for f in `ls Gene_*txt_unique`
do
	chr_num=`awk '{print $2}' ${f} | cut -d'_' -f 1 | head -n1`
	outfn=SNP_loc_${chr_num}_100bp_middle.bed
	paste <(cat ${f} | awk '{print $2}' | awk -F'_' '{print $1}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=50 '{print $2-s}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=50 '{print $2+s}')  <(cat ${f} | awk '{print $2}') >> ${outputdir}/${outfn}
        outfn=SNP_loc_${chr_num}_100bp_1kb_left.bed
        paste <(cat ${f} | awk '{print $2}' | awk -F'_' '{print $1}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=500 '{print $2-s}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=50 '{print $2-s}')  <(cat ${f} | awk '{print $2}') >> ${outputdir}/${outfn}
        outfn=SNP_loc_${chr_num}_100bp_1kb_right.bed
        paste <(cat ${f} | awk '{print $2}' | awk -F'_' '{print $1}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=50 '{print $2+s}') <(cat ${f} | awk '{print $2}' | awk -F'_' -v s=500 '{print $2+s}')  <(cat ${f} | awk '{print $2}') >> ${outputdir}/${outfn}
done



cd ${outputdir}
fn=commands.txt
for f in `ls SNP_loc_chr*bed`
do
	echo "sort -k1,1 -k2,2n ${f} | uniq > temp_${f}" >> ${fn}
done

ml parallel
parallel :::: ${fn}


for f in `ls temp_*`
do
	mv ${f} ${f#temp_}
done	

