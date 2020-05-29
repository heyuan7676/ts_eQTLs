#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00
#SBATCH -p debug

alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
#alignmentFn=${alignmetn_dir}/liver_FOXA1_ChIP_seq_noWASPAligned.sortedByCoord.out.bam
alignmentFn=${alignmetn_dir}/liver_FOXA1_ChIP_seqAligned.sortedByCoord.out.bam


### reads that overlap with variants & pass WASP filter
ml samtools
samtools view ${alignmentFn} | grep vA | grep vW:i:1 | grep -v "vA:B:c,3" > ${alignmetn_dir}/tp.sam


# run second chunk of 5_TFBS_ChIP_ASB.ipynb

## combine
#outfn=${alignmentFn}.variants.bed
#paste <(cat ${alignmetn_dir}/bedfile_col1_3.bed) <(cat ${alignmetn_dir}/bedfile_col1_4.bed) | sed 's/ /	/g' > ${outfn}
#awk '{print $1"_"$3, $4,$1}' ${outfn} | sed 's/ /	/g' > ${outfn}_SNPname

#rm ${alignmetn_dir}/tp*
#rm ${alignmetn_dir}/bedfile_col*bed
