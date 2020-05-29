#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --no-requeue
#SBATCH -p lrgmem

#bash vcfFn_process.sh

## compute LD blocks
module load plink
r2Thr="$1"
chr="$2"

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes
vcffn=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01
outdir=${datadir}/plink_r2_${r2Thr}_bychr

mkdir -p ${outdir}
mkdir -p ${datadir}/plink_chr

chr=chr${chr}
echo ${chr}

## transform VCF files to plink binary files
vcfFn_chr=${datadir}/vcf_chr/${vcffn}_${chr}.vcf
binaryFn=${datadir}/plink_chr/genotypes_v8_${chr}
plink --vcf ${vcfFn_chr} --recode --out ${binaryFn}

## deal with the problem of SNP name!!!!
paste <(awk '{print $1}' ${binaryFn}.map) <(awk '{print $2}' ${binaryFn}.map| cut -d'_' -f1,2 ) <(awk '{print $3,$4}' ${binaryFn}.map) > ${binaryFn}.map.new
sed 's/ /	/g' ${binaryFn}.map.new > ${binaryFn}.map
rm ${binaryFn}.map.new

## compute R2
r2Out=${outdir}/plink_r2_${r2Thr}_${chr}
plink --file ${binaryFn} --extract ${datadir}/SNPIDs.txt --r2 inter-chr --ld-window-r2 ${r2Thr} --out ${r2Out}

