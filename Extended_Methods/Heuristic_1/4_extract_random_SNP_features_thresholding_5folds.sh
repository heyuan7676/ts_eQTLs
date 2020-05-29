#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --no-requeue

### produce:  SNPs x Tissues in DNase data


ml parallel

source ./GLOBAL_VAR.sh


iddir=${pairdir}/match_random/batches_random/
command_fn=extractSNP_features_commands_closeToTop_5folds.txt
rm -f ${command_fn}

echo "extract SNP features"
for tis in `cat tissues.txt`
do
                cd ${iddir}
		prefix=${tis}_ts_ciseQTL_closeToTop
                fns=`ls ${prefix}*idx*.txt  | head -n5`
                idfn=${prefix}_outlierPairs_random_matched_5folds.txt
                rm -f ${idfn}
                cat ${fns} >> ${idfn}

		cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/3_analysis
		for anno_feature in DNase ROADMAP
		do
			featurefn=SNP_loc_${anno_feature}.bed
			outfn=Thresholding_outlierSNPs_random_matched_${anno_feature}_${tis}_5folds.txt
			echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
		done

done


tis=Shared
                cd ${iddir}
                prefix=${tis}_ts_ciseQTL_closeToTop
                fns=`ls ${prefix}*idx*.txt  | head -n5`
                idfn=${prefix}_outlierPairs_random_matched_5folds.txt
                rm -f ${idfn}
                cat ${fns} >> ${idfn}

                cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/3_analysis
                for anno_feature in DNase ROADMAP
                do
                        #echo ${anno_feature}
                        featurefn=SNP_loc_${anno_feature}.bed
			outfn=Thresholding_outlierSNPs_random_matched_${anno_feature}_${tis}_5folds.txt
                        echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
                done



parallel :::: ${command_fn}
