#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -p quickdebug

featureName=ENCODE
#featureName=DNase
datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/prediction_model/tissue_${featureName}
cd ${datadir}


#for subfeature in  H3K4me1 H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K9me3
for subfeature in H3K36me3 H3K9me3
#for subfeature in DNase
do
	echo $subfeature
	for loc in 100bp_1kb_left 100bp_1kb_right 100bp_middle
	do
        	echo $loc
		fn=SNP_loc_${loc}_${subfeature}.bed
		newfn=SNP_loc_${loc}_${subfeature}_onlyhits.bed
		cat ${fn} | cut -d'	' -f2-  | sed 's/,//g' | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}' > temp

		paste temp ${fn} > ${fn}_temp
		head -n1 ${fn} > ${newfn}
		awk '{if($1!=0) print $0}' ${fn}_temp | cut -d'	' -f2- >> ${newfn}

		fn=SNP_loc_${loc}_${subfeature}_coverage.bed
		newfn=SNP_loc_${loc}_${subfeature}_coverage_onlyhits.bed
                paste temp ${fn} > ${fn}_temp
                head -n1 ${fn} > ${newfn}
                awk '{if($1!=0) print $0}' ${fn}_temp | cut -d'	' -f2- >> ${newfn}

	done


	for loc in TES TSS_500bp_1kb TSS_500bp_left TSS_500bp_right
        do
                echo $loc
                fn=gene_loc_${loc}_${subfeature}.bed
                newfn=gene_loc_${loc}_${subfeature}_onlyhits.bed
                cat ${fn} | cut -d'	' -f2-  | sed 's/,//g' | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}' > temp
	
                paste temp ${fn} > ${fn}_temp
                head -n1 ${fn} > ${newfn}
                awk '{if($1!=0) print $0}' ${fn}_temp | cut -d'	' -f2- >> ${newfn}

                fn=gene_loc_${loc}_${subfeature}_coverage.bed
                newfn=gene_loc_${loc}_${subfeature}_coverage_onlyhits.bed
                paste temp ${fn} > ${fn}_temp
                head -n1 ${fn} > ${newfn}
		awk '{if($1!=0) print $0}' ${fn}_temp | cut -d'	' -f2- >> ${newfn}
 
        done
done

rm *_temp
