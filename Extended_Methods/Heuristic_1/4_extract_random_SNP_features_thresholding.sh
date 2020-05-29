#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue
#SBATCH -p shared
 
### produce:  SNPs x Tissues in DNase data


ml parallel

source ./GLOBAL_VAR.sh


iddir=${pairdir}/match_random/batches_random/
command_fn=extractSNP_features_commands_closeToTop.txt
rm -f ${command_fn}

echo "extract SNP features"
#for tis in `cat tissues.txt`
#do
tis=Shared
                for idx in {0..100}
                do
                        idfn=${tis}_ts_ciseQTL_closeToTop_random_matched_idx${idx}.txt
                        for anno_feature in DNase ROADMAP
                        do
                                featurefn=SNP_loc_${anno_feature}.bed
                                outfn=Thresholding_outlierSNPs_random_matched_${anno_feature}_${tis}_idx${idx}.txt
                                rm -f ${sigSNPfeaturedir}/${outfn}

				#if [[ -f "${sigSNPfeaturedir}/${outfn}" ]]
                                #then
					#L1=`awk '{print $2}' ${iddir}/${idfn} | sort | head -n1`
                                        #L2=`sed "1d" ${sigSNPfeaturedir}/${outfn} | awk '{print $1}' | sort | head -n1`
                                        #L1=`awk '{print $2}' ${iddir}/${idfn} | sort | uniq | wc -l | awk '{print $1}'`
                                        #L2=`sed "1d" ${sigSNPfeaturedir}/${outfn} | wc -l | awk '{print $1}'`
					#if [ "$L1" -eq "$L2" ]
                                        #then
                                        #        continue
                                        #fi
                                        # exist, but not finished
                                        #echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
                                        #continue
                                #fi
				echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
			done
		done

#done



parallel :::: ${command_fn}
