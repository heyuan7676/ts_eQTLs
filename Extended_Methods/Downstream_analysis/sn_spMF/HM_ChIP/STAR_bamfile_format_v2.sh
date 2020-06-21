#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00
#SBATCH -p debug

ml samtools

for hm in H3K4me1 H3K4me3 H3K9me3 H3K27ac H3K27me3 H3K36me3
do
	echo $hm

	alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/HM_ChIP_seq/STAR_output
	alignmentFn=${alignmetn_dir}/liver_${hm}_seqAligned.sortedByCoord.out.bam

	samtools view ${alignmentFn} | grep vA | grep vW:i:1 | grep -v "vA:B:c,3" > ${alignmetn_dir}/${hm}.sam.txt
done

