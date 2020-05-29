#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p lrgmem


### get the intersection of SNPs in TFBS and SNPs in the outlier Pairs
alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
alignmentFn=liver_HNF4A_ChIP_seqAligned.sortedByCoord.out.bam
SNP_in_TFBS=${alignmetn_dir}/${alignmentFn}.variants.bed_SNPname


output_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_outlierSNPs
cd ${output_dir}

### deal with chromsome first!!!! -- because the vcf files are too large, can't index each SNP one by one
rm chr*
awk '{print >> $3.txt}' ${SNP_in_TFBS}


### get the SNP info from VCF file
SNP_vcf_chr_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
output_fn_merged=SNP_inTFBS_inOutlier_REF_AL.txt
rm -f ${output_fn_merged}
for chromosome in {1..22}
do
	echo chromosome${chromosome}
	cat chr${chromosome} | uniq | sort -k2,2 | sed  's/ /	/g'  > temp
	mv temp chr${chromosome} 
	vcfn=${SNP_vcf_chr_dir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chromosome}.vcf
	cat ${vcfn} | grep -v "#" | cut -f3-5 | awk '{split($0,a,"_");print a[1]"_"a[2],$0}'  | awk 'NR==FNR {id[$1]; next} $1 in id'  chr${chromosome} - | sort -k1,1 | sed 's/ /	/g' > chr${chromosome}.tp
	join -1 1 -2 2 chr${chromosome}.tp chr${chromosome} >> ${output_fn_merged}
done

#rm chr*
#rm tp

### map back to groups
#pairdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets/
#LMfn=SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110_LD1_Loadings_beta_BH_corrected_alpha0.05
#for group in {0..22}
#do
#	echo group${group}
#        outfn_final=${output_dir}/${LMfn}_group${group}_${alignmentFn%.bam}.txt
#	echo "SNP_name	SNPID_vcf	REF	AL	TFBS_SNP"	> ${outfn_final}
#	awk 'NR==FNR {id[$2]; next} $1 in id'  ${pairdir}/${LMfn}_outlierPairs${randomString}_group${group}.txt ${output_fn_merged} | uniq | sed 's/ /	/g' >> ${outfn_final}
#done







