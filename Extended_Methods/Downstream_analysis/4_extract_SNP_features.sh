#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --no-requeue
#SBATCH -p lrgmem

### produce:  SNPs x Tissues in DNase data

source ./GLOBAL_VAR.sh

function extractSNP_feature(){
        iddir=$1
        idfn=$2
        inputdir=$3
        inputfn=$4

        outdir=$5
        outfn=$6

        head -n1 ${inputdir}/${inputfn} > ${outdir}/${outfn}
        awk -F '\t' 'NR==FNR {id[$2]; next} $1 in id' ${iddir}/${idfn} ${inputdir}/${inputfn} >> ${outdir}/${outfn}

        cat ${outdir}/${outfn} | sed "s/\.\t/0\t/g" | sed "s/\.$/0/" | sed 's/ /	/g'  > ${outdir}/${outfn}_temp
        mv ${outdir}/${outfn}_temp ${outdir}/${outfn}
        echo ${outdir}/${outfn}
        awk '{print NF}' ${outdir}/${outfn} | uniq
}



iddir=${pairdir}/match_random
dir=`pwd`

echo "extract SNP features"
for group in {0..22}
do
	echo group$group

	cd ${iddir}/batches_random/
	fns=`ls ${LMfn}*group${group}_*idx*.txt  | head -n4`
	idfn_rd=${LMfn}_outlierPairs_random_matched_group${group}_4folds.txt
	rm -f ${iddir}/${idfn_rd}
	cat ${fns} >> ${iddir}/${idfn_rd}

	cd ${dir}
	for anno_feature in DNase ROADMAP
	do
		featurefn=SNP_loc_${anno_feature}.bed

		outfn=${FMfn}_outlierSNPs_${anno_feature}_group${group}.txt
		idfn=${LMfn}_outlierPairs_group${group}_annotated.txt
		extractSNP_feature ${iddir} ${idfn} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}

		outfn=${FMfn}_outlierSNPs_random_matched_${anno_feature}_group${group}_4folds.txt
		idfn_rd=${LMfn}_outlierPairs_random_matched_group${group}_4folds.txt
		extractSNP_feature ${iddir} ${idfn_rd} ${allSNPfeaturedir}/tissue_${anno_feature} ${featurefn} ${sigSNPfeaturedir} ${outfn}
	done

done


