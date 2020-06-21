#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

TF="$1"

fn=ASB_SNPs_${TF}_Filter10_FDR005.txt 
out_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/ChIP_ASB/Direction


#### VCF file for SNPs
vcf_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
SNP_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/ChIP_ASB/ASB_SNPS

outfn=${fn%.txt}_SNPs_vcf.txt
vcf_fn=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr1.vcf
cat ${vcf_dir}/${vcf_fn} | grep -v "##" | grep "GTEX-" > ${out_dir}/${outfn}

for chromosome in {1..22}
do
	echo $chromosome
	vcf_fn=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chromosome}.vcf
	grep -Ff ${SNP_dir}/${fn} ${vcf_dir}/${vcf_fn} >> ${out_dir}/${outfn}
done



#### TPM file for eGenes
eQTL_fn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets/SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110_LD1_Loadings_beta_BH_corrected_alpha0.05_outlierPairs_group15.txt
TPM_fn=/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/Liver.txt

grep -Ff ${SNP_dir}/${fn} ${eQTL_fn} | awk '{print $1}' > temp_${fn}
head -n1 $TPM_fn > ${out_dir}/${fn%.txt}_genes_TPM.txt
grep -Ff temp_${fn} $TPM_fn >> ${out_dir}/${fn%.txt}_genes_TPM.txt
rm temp_${fn}




#### ChIP-seq reads for ASB
fastq_reads_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
for chipseq_fn in `ls ${fastq_reads_dir}/*${TF}*formatted | grep -v noWASP | xargs -n1 basename`
do
	echo $chipseq_fn
	grep -Ff ${SNP_dir}/${fn} ${fastq_reads_dir}/${chipseq_fn} > ${out_dir}/${chipseq_fn%_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_formatted}.txt
done





