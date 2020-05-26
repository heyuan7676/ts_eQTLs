#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

source ./GLOBAL_VAR.sh

featureName=ROADMAP_"$1"
anno_feature=TFmotif

function extractSNP_feature(){
        iddir=$1
        idfn=$2
        inputdir=$3
        inputfn=$4

        outdir=$5
        outfn=$6

        head -n1 ${inputdir}/${inputfn} > ${outdir}/${outfn}
        awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' ${iddir}/${idfn} ${inputdir}/${inputfn} >> ${outdir}/${outfn}

        cat ${outdir}/${outfn} | sed "s/\.\t/0\t/g" | sed "s/\.$/0/" > ${outdir}/${outfn}_temp
        mv ${outdir}/${outfn}_temp ${outdir}/${outfn}
}



for group in {0..22}
do
	echo $group
	group=group${group}
	iddir=${activeSNPdir}
	idfn=${FMfn}_Active_SNPset_${featureName}_${group}.txt

	featurefn=SNP_loc_${anno_feature}.bed
	outfn=${FMfn}_Active_SNPset_${featureName}_${anno_feature}_${group}.txt
	extractSNP_feature ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${activeSNPfeaturedir} ${outfn}

done

