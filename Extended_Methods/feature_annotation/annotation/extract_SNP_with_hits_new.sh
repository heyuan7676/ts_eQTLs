#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

subfeature="$1"
	echo $subfeature
	cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/allPairs/tissue_${subfeature}
	fn=SNP_loc_${subfeature}.bed
	newfn=SNP_loc_${subfeature}_onlyhits.bed
	
	cat ${fn} | cut -d'	' -f2-  | sed 's/,//g' | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}' > temp

	paste temp ${fn} > ${fn}_temp
	head -n1 ${fn} > ${newfn}
	awk '{if($1!=0) print $0}' ${fn}_temp | cut -d'	' -f2- >> ${newfn}

