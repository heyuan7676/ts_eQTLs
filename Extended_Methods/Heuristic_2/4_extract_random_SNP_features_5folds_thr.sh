#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue
#SBATCH -p shared

### produce:  SNPs x Tissues in DNase data


ml parallel

source ./GLOBAL_VAR.sh
command_fn=command_spMF.txt
rm -f ${command_fn}


iddir=${pairdir}/match_random/batches_random/
prefix=${LMfn}
echo "extract SNP features"
for group in {-1..27}
do
		cd ${iddir}
		fns=`ls ${prefix}*group${group}_*idx*.txt  | head -n5`
		idfn=${prefix}_outlierPairs_random_matched_group${group}_5folds.txt
		rm -f ${idfn}
		cat ${fns} >> ${idfn}
		cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/revision/Thresholding
		for anno_feature in DNase ROADMAP 
		do
			featurefn=SNP_loc_${anno_feature}.bed
			outfn=${prefix}_outlierSNPs_random_matched_${anno_feature}_group${group}_5folds.txt
			rm -f ${outfn}
			echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
		done

done


parallel :::: ${command_fn}
