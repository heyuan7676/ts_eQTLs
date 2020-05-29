#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=2:00:00
#SBATCH -p quickdebug

source ./GLOBAL_VAR.sh

fn=annotated_commands.txt
rm -f ${fn}
#featureName=ENCODE
featureName=DNase
#for subfeature in  H3K4me1 H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K9me3
for subfeature in DNase
do
echo $subfeature
featuredir=${datadir}/annotations/${featureName}/merged
#cd ${featuredir}
#for featurefn in `ls *${subfeature}*bed`
#do
        #echo $featurefn
        #sort -k1,1 -k2,2n $featurefn > tmp_$featurefn
        #mv tmp_$featurefn $featurefn
#done



cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/0_feature_annotation/annotation
for loc in 100bp_1kb_left 100bp_1kb_right 100bp_middle
do
	echo $loc
	for chr in {1..22}
	do
		echo "bash annotate_SNP_gene.sh ${featureName} ${subfeature} SNP_loc_chr${chr}_${loc}.bed SNP_loc" >> ${fn}
	done
done


for loc in TES TSS_500bp_1kb TSS_500bp_left TSS_500bp_right
do
	echo $loc
	echo "bash annotate_SNP_gene.sh ${featureName} ${subfeature} gene_loc_${loc}.bed gene_loc" >> ${fn}
done
 

ml parallel
parallel :::: ${fn}



cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/prediction_model/tissue_${featureName}
for loc in 100bp_1kb_left 100bp_1kb_right 100bp_middle
do
	fn=SNP_loc_${loc}_${subfeature}.bed
	cat SNP_loc_chr1_${loc}_${subfeature}.bed >  ${fn}
        for chr in {2..22}
        do
		sed "1d" SNP_loc_chr${chr}_${loc}_${subfeature}.bed >> ${fn}
        done

	fn=SNP_loc_${loc}_${subfeature}_coverage.bed
        cat SNP_loc_chr1_${loc}_${subfeature}_coverage.bed >  ${fn}
        for chr in {2..22}
        do
                sed "1d" SNP_loc_chr${chr}_${loc}_${subfeature}_coverage.bed >> ${fn}
        done

done
done
