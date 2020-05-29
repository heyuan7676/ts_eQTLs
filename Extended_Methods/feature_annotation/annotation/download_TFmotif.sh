#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif
wget -r --no-parent http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/tsv/
