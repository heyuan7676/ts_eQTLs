#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --no-requeue
#SBATCH -p debug

cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
fn=maf_commands.txt
rm -r ${fn}
for chr in {1..22}
do
	echo "cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chr}.vcf | grep ^[^#] | cut -d'	' -f1,2 | sed 's/	/_/g' > chr${chr}_ID.txt" >> ${fn}
	echo "cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_chr${chr}.vcf | grep ^[^#] | cut -d'	' -f8 | cut -d';' -f2  | cut -d'=' -f2 > chr${chr}_AF.txt" >> ${fn}
done

ml parallel
parallel :::: ${fn}


for i in {1..22}
do
	paste <(cat chr${i}_ID.txt) <(cat chr${i}_AF.txt) > chr${i}_SNP_AF.txt 
done


cat chr*_SNP_AF.txt >> v8_SNP_AF.txt
rm chr*
