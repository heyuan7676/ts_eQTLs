#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p shared

ml MACS2

peaks_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks
data_dir=${peaks_dir}/output

fn=CTCF_S1


cd /home-4/yhe23@jhu.edu/work/yuan/tools/UCSC_tools 

bash bdg2bw.sh ${data_dir}/${fn}_FE.bdg hg38.len 
#bdg2bw ${data_dir}/${fn}_logLR.bdg hg38.len

