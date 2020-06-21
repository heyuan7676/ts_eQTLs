#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p shared

ml MACS2
ml samtools

peaks_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks

fn="$1"
macs2 callpeak -t ${peaks_dir}/input/${fn}.bam -c ${peaks_dir}/input/${fn}_control.bam -B -n ${fn}_peaks --outdir ${peaks_dir}/output -f BAM -g hs -q 0.1



