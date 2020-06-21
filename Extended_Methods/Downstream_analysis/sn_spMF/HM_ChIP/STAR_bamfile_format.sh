#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00
#SBATCH -p debug

alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
#alignmentFn=${alignmetn_dir}/liver_HNF4A_ChIP_seq_noWASPAligned.sortedByCoord.out.bam
alignmentFn=${alignmetn_dir}/liver_HNF4A_ChIP_seqAligned.sortedByCoord.out.bam


### reads that overlap with variants & pass WASP filter
ml samtools
samtools view ${alignmentFn} | grep vA | grep vW:i:1 | grep -v "vA:B:c,3" > ${alignmetn_dir}/tp.sam


### extract the peaks mapped locations
paste <(cat ${alignmetn_dir}/tp.sam  | awk '{print $3}') <(cat ${alignmetn_dir}/tp.sam  | awk '{print $19}' | cut -d"," -f2) > ${alignmetn_dir}/tp.sam1
awk -v s=1 '{print "chr"$1, $2, $2+1}' ${alignmetn_dir}/tp.sam1 > ${alignmetn_dir}/bedfile_col1_3.bed

### get the variant
paste <(cat ${alignmetn_dir}/tp.sam  | awk '{print $4}') <(cat ${alignmetn_dir}/tp.sam  | awk '{print $19}' | cut -d"," -f2) > ${alignmetn_dir}/tp.sam3
paste <(awk '{print $10}' ${alignmetn_dir}/tp.sam) <(awk '{ $3 = $2 - $1 } 1' ${alignmetn_dir}/tp.sam3 | awk '{print $3}') > ${alignmetn_dir}/tp.sam4

# run first chunk of 5_TFBS_ChIP_ASB.ipynb
#variant_fn=${alignmetn_dir}/bedfile_col1_4.bed
#rm -f ${variant_fn}
#while IFS= read -r line
#do 
#	array=`echo $line | awk '{print $1}'`
#	loc=`echo $line | awk '{print $2}'`
#	loc=$(( loc + 1 ))
#	echo ${array} | cut -c${loc} >> ${variant_fn}
#done < ${alignmetn_dir}/tp.sam4


## combine
outfn=${alignmentFn}.variants.bed
paste <(cat ${alignmetn_dir}/bedfile_col1_3.bed) <(cat ${alignmetn_dir}/bedfile_col1_4.bed) | sed 's/ /	/g' > ${outfn}
awk '{print $1"_"$3, $4,$1}' ${outfn} | sed 's/ /	/g' > ${outfn}_SNPname

#rm ${alignmetn_dir}/tp*
#rm ${alignmetn_dir}/bedfile_col*bed
