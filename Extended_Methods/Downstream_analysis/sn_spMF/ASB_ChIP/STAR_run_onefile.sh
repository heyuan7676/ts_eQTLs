#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00
#SBATCH -p lrgmem

## fn: the file that has the fastq file in /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/fastq/

#fn=remap.HNF4A.liver.fastq.file.txt
#fn="$1"
#bash STAR_AL_alignment.sh ${fn}

#fn=${fn%.fastq.file.txt}_ChIP_seqAligned.sortedByCoord.out.bam

fn="$1"
bash STAR_bamfile_format.sh ${fn}


fn=${fn}.filtered_noWASP.txt
bash STAR_intersect_outlier_SNPs.sh ${fn}

