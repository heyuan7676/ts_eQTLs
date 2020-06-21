#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue
#SBATCH -p shared
 
### produce:  SNPs x Tissues in DNase data


ml parallel

source ./GLOBAL_VAR.sh


iddir=${pairdir}
command_fn=extractSNP_features_commands_closeToTop.txt
rm -f ${command_fn}

echo "extract SNP features"
for tis in `cat tissues.txt`
do
		idfn=${tis}_ts_ciseQTL_closeToTop.txt
		for anno_feature in DNase ROADMAP
		do
			#echo ${anno_feature}
			featurefn=SNP_loc_${anno_feature}.bed
			outfn=Thresholding_outlierSNPs_${anno_feature}_${tis}.txt
			echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
		done

done


tis=Shared
                idfn=${tis}_ts_ciseQTL_closeToTop.txt
                for anno_feature in DNase ROADMAP
                do      
                        #echo ${anno_feature}
                        featurefn=SNP_loc_${anno_feature}.bed
			outfn=Thresholding_outlierSNPs_${anno_feature}_${tis}.txt
                        echo "bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}" >> ${command_fn}
                done


parallel :::: $command_fn