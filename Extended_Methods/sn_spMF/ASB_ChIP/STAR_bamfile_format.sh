#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=20:00:00
#SBATCH -p lrgmem


ml samtools
ml parallel

alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
cd ${alignmetn_dir}
fn="$1"



#rm -f command_fn.txt
#for fn in `ls *bam`
#do
	### reads that overlap with variants & pass WASP filter
	#samtools view ${fn} | grep vA | grep vW:i:1 | grep -v "vA:B:c,3" > ${fn}.filtered.txt
	samtools view ${fn} | grep vA  > ${fn}.filtered_noWASP.txt
#done
#parallel :::: command_fn.txt

