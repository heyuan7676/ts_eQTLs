#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p express

ml MACS2

peaks_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks
data_dir=${peaks_dir}/output

fn="$1"
macs2 bdgcmp -t ${data_dir}/${fn}_peaks_treat_pileup.bdg -c ${data_dir}/${fn}_peaks_control_lambda.bdg -o ${data_dir}/${fn}_FE.bdg -m FE
macs2 bdgcmp -t ${data_dir}/${fn}_peaks_treat_pileup.bdg -c ${data_dir}/${fn}_peaks_control_lambda.bdg -o ${data_dir}/${fn}_logLR.bdg -m logLR -p 0.00001



cd /home-4/yhe23@jhu.edu/work/yuan/tools/UCSC_tools

bash bdg2bw.sh ${data_dir}/${fn}_FE.bdg hg38.len 
bash bdg2bw.sh ${data_dir}/${fn}_logLR.bdg hg38.len

