#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --no-requeue
#SBATCH -p quickdebug

module load plink
module load htslib

vcfdir=/work-zfs/abattle4/lab_data/GTEx_v8/genotypes_5_22
vcffn=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01
inputFn=${vcfdir}/${vcffn}.vcf.gz

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes
r2Thr=0.2

for chr in {1..22}
do
	chr=chr${chr}
	echo $chr

	## split VCF file
	vcfFn_chr=${datadir}/vcf_chr/${vcffn}_${chr}.vcf.gz
	tabix -h ${inputFn} ${chr} > ${vcfFn_chr}
done

