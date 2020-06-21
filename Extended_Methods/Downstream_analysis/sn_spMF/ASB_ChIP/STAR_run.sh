#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

#for TF in HNF4A CTCF MAX JUND ZBTB33
#do
#	for fn in `ls /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/fastq/${TF}_Liver*.fastq.file.txt`
#	do
#		fn=`echo ${fn} | xargs -n1 basename`
		#echo $fn
		#sbatch STAR_AL_alignment.sh ${fn}
	#done
#done


#bash STAR_bamfile_format.sh



#echo "Step 2) .formatted" and 3) GTEx SNPs (SNPname)"
#alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
#for fn in `ls ${alignmetn_dir}/*.filtered.txt | xargs -n1 basename`
#do
#        sbatch STAR_intersect_outlier_SNPs.sh ${fn}
#done	


alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
for alignmentFn in `ls ${alignmetn_dir}/MAX*formatted | xargs -n1 basename `
do
	echo $alignmentFn
	bash ASB_run.sh ${alignmentFn}
done
