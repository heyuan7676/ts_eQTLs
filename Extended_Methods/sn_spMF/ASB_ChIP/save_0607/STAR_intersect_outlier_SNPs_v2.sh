#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


### get the intersection of SNPs in TFBS and SNPs in the outlier Pairs
alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
alignmentFn=liver_FOXA1_ChIP_seqAligned.sortedByCoord.out.bam
SNPname=${alignmetn_dir}/${alignmentFn}.formatted


output_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_outlierSNPs_v2
cd ${output_dir}

### deal with chromsome first!!!! -- because the vcf files are too large, can't index each SNP one by one
rm -f chr*
awk '{print >> "chr"$1.txt}' ${SNPname}


### get the SNP info from VCF file
SNP_vcf_chr_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
output_fn_merged=SNP_inTFBS_inOutlier_REF_AL_FOXA1.txt
rm -f ${output_fn_merged}
for chromosome in {1..22}
do
	echo chromosome${chromosome}
	cat chr${chromosome} | uniq | sort -k4,4 | sed  's/ /	/g'  > temp
	mv temp chr${chromosome} 
	vcfn=${SNP_vcf_chr_dir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chromosome}.vcf
	cat ${vcfn} | grep -v "#" | cut -f3-5 | awk '{split($0,a,"_");print a[1]"_"a[2],$0}'  | awk 'NR==FNR {id[$4]; next} $1 in id'  chr${chromosome} - | sort -k1,1 | sed 's/ /	/g' > chr${chromosome}.tp
	join -1 1 -2 4 chr${chromosome}.tp chr${chromosome} >> ${output_fn_merged}
done








