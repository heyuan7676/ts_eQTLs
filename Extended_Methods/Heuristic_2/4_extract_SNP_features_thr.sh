#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --no-requeue
#SBATCH -p lrgmem

### produce:  SNPs x Tissues in DNase data


ml parallel

source ./GLOBAL_VAR.sh
rm -f ${command_fn}


iddir=${pairdir}/match_random
command_fn=extractSNP_features_commands.txt

echo "extract SNP features"
for group in {-1..27}
do
	echo group$group
		idfn=${LMfn}_outlierPairs_group${group}_annotated.txt
		for anno_feature in DNase ROADMAP
		do
			#echo ${anno_feature}
			featurefn=SNP_loc_${anno_feature}.bed
			outfn=${LMfn}_outlierSNPs_${anno_feature}_group${group}.txt
			rm -f ${outfn}
			bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}
		done

done


