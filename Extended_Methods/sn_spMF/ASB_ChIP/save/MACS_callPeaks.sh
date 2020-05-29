#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH -p shared

ml MACS2

bamfile_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
#fn=HNF4A_Liver_single-ended_ENCFF302XOK_ChIP_seqAligned.sortedByCoord.out.bam
fn="$1"

macs2 callpeak -t ${bamfile_dir}/${fn} -n ${fn%_ChIP_seqAligned.sortedByCoord.out.bam} --outdir /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks -f BAM -g hs -q 0.1
