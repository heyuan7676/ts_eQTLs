#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

source ./GLOBAL_VAR.sh
ml parallel

featureName=ROADMAP_"$1"
anno_feature=TFmotif
fn=extract_active_SNP_features_commands_${featureName}.txt
rm -f ${fn}


rm -f ${activeSNPfeaturedir}/*${featureName}_${anno_feature}*
for group in {-1..27}
do
	echo $group
	group=group${group}
	iddir=${activeSNPdir}
	idfn=${FMfn}_Active_SNPset_${featureName}_${group}.txt

	featurefn=SNP_loc_${anno_feature}.bed
	outfn=${FMfn}_Active_SNPset_${featureName}_${anno_feature}_${group}.txt
	bash 5_extractSNP_feature_2.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${activeSNPfeaturedir} ${outfn}
	#echo "bash 5_extractSNP_feature_2.sh ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${activeSNPfeaturedir} ${outfn}" >> ${fn}

done
#parallel :::: ${fn}

