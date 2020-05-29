#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH -p debug


function indexPairsFromciseQTL {
	tissue="$1"
        cisFn=${cisDir}/${tissue}.allpairs.txt
	inputFn=${cisFn_matched_SNPs}/${tissue}.allpairs.txt
        paste <(cat ${cisFn} | cut -d'_' -f 1,2) <(cat ${cisFn} | cut -d'	' -f3-) > ${inputFn}

	outFn=${outDir}/${tissue}.v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.txt
	rm -f ${outFn}
	awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' ${idsFn} ${inputFn} | sort -u  | sed 's/ /	/g' | awk '!seen[$1$2]++' > ${outFn}

	N=`wc -l ${outFn} | awk '{print $1}'`
	echo 'Tissue:' $tissue '-' ${N}

	leftoutFn=${outDir}/${tissue}.leftout.txt
	rm -f ${leftoutFn}
	comm -13 <(awk '{print $1,$2}' ${outFn} | sed 's/ /	/g') <(sed "1d" ${idsFn} | awk '{print $1,$2}' | sed 's/ /	/g') > ${leftoutFn}
	wc -l ${leftoutFn}

	completeFn=${outDir}/${tissue}.${suffix}.txt
	cat ${outFn}  > ${outDir}/temp_${tissue}
	cat ${leftoutFn} >> ${outDir}/temp_${tissue}
	
	head -n1 ${inputFn} > ${completeFn}
	cat ${outDir}/temp_${tissue} | sort >> ${completeFn}
	rm -f temp_${tissue}
	wc -l ${completeFn}
}

