#!/bin/bash

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




extractSNP_feature "$1" "$2" "$3" "$4" "$5" "$6"
