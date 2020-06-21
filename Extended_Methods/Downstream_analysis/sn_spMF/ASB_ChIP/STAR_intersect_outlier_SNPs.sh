#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


fn="$1"

ml python/2.7
python 5_TFBS_ChIP_ASB_format_samFiles.py ${fn}


### get the intersection of SNPs in TFBS and SNPs in the outlier Pairs
alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
output_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output_GTExSNPs

alignmentFn=${fn}_formatted

rm -f -r ${alignmetn_dir}/${alignmentFn}_folder
mkdir ${alignmetn_dir}/${alignmentFn}_folder
mv ${alignmetn_dir}/${alignmentFn} ${alignmetn_dir}/${alignmentFn}_folder/

cd ${alignmetn_dir}/${alignmentFn}_folder/

### deal with chromsome first!!!! -- because the vcf files are too large, can't index each SNP one by one
rm -f chr*
awk '{print >> "chr"$2.txt}' ${alignmentFn}

### get the SNP info from VCF file
SNP_vcf_chr_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
output_fn_merged=${output_dir}/SNP_inTFBS_inGTEx_${alignmentFn}
rm -f ${output_fn_merged}
for chromosome in {1..22}
do
	echo chromosome${chromosome}
	cat chr${chromosome} | uniq | sort -k1,1 | sed  's/ /	/g'  > temp
	mv temp chr${chromosome} 
	vcfn=${SNP_vcf_chr_dir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chromosome}.vcf
	cat ${vcfn} | grep -v "#" | cut -f3-5 | awk '{split($0,a,"_");print a[1]"_"a[2],$0}'  | awk 'NR==FNR {id[$1]; next} $1 in id'  chr${chromosome} - | sort -k1,1 | sed 's/ /	/g' > chr${chromosome}.tp
	join -1 1 -2 1 chr${chromosome}.tp chr${chromosome} >> ${output_fn_merged}
done

mv ${alignmentFn} ${alignmetn_dir}/
rm -r ${alignmetn_dir}/${alignmentFn}_folder




