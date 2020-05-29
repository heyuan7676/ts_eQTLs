#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=2:00:00
#SBATCH -p quickdebug

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets

inputfn=/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt
outputdir=${datadir}/prediction_model/gene_loc



rm -r -f ${outputdir}
mkdir -p ${outputdir}


outfn=gene_loc_TSS_500bp.bed
paste <(sed 1d ${inputfn} | awk '{print $2}') <(sed 1d ${inputfn} | awk -v s=500 '{print $3-s}') <(sed 1d ${inputfn} | awk '{print $3}')  <(sed 1d ${inputfn} | awk '{print $1}') | grep -v "chrM" | grep -v "chrX" | grep -v "chrY" >> ${outputdir}/${outfn}

outfn=gene_loc_TSS_500bp_1kb.bed
paste <(sed 1d ${inputfn} | awk '{print $2}') <(sed 1d ${inputfn} | awk -v s=1000 '{print $3-s}') <(sed 1d ${inputfn} | awk -v s=500 '{print $3-s}')  <(sed 1d ${inputfn} | awk '{print $1}')  | grep -v "chrM" | grep -v "chrX" | grep -v "chrY">> ${outputdir}/${outfn}

outfn=gene_loc_TSS_500bp_right.bed
paste <(sed 1d ${inputfn} | awk '{print $2}') <(sed 1d ${inputfn} | awk '{print $3}') <(sed 1d ${inputfn} | awk -v s=500 '{print $3+s}')  <(sed 1d ${inputfn} | awk '{print $1}')  | grep -v "chrM" | grep -v "chrX" | grep -v "chrY" >> ${outputdir}/${outfn}

outfn=gene_loc_TES.bed
paste <(sed 1d ${inputfn} | awk '{print $2}') <(sed 1d ${inputfn} | awk '{print $4}') <(sed 1d ${inputfn} | awk -v s=1000 '{print $4+s}')  <(sed 1d ${inputfn} | awk '{print $1}')  | grep -v "chrM" | grep -v "chrX" | grep -v "chrY" >> ${outputdir}/${outfn}


cd ${outputdir}
fn=commands.txt
for f in `ls gene_loc*bed`
do
	echo "sort -k1,1 -k2,2n ${f} | uniq > temp_${f}" >> ${fn}
done

ml parallel
parallel :::: ${fn}


for f in `ls temp_*`
do
	mv ${f} ${f#temp_}
done	

