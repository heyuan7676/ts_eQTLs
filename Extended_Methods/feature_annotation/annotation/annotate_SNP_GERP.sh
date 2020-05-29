#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00

source ./GLOBAL_VAR.sh

feature=GERP
featuredir=${datadir}/annotations/GERP

inputdir=${datadir}/annotations/cbset_pairs
inputfn=${suffix}.SNP_loc.txt
outputdir=${inputdir}/tissue_${feature}

mkdir -p ${outputdir}

function rename_uniqueRows(){
        fn="$1"
        y=${fn/(var.2)/-var.2}
        y=${y/(var.3)/-var.3}
        cat ${fn} |  awk '!seen[$4]++' | awk '{print $4, $8} ' > ${y}_unique
	mv ${y}_unique ${y}
}

featurefn=GRCh38_GERP.bed
outFn=${inputfn%.txt}_${featurefn}
/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wao -a ${inputdir}/${inputfn} -b ${featuredir}/${featurefn} -nonamecheck -sorted | sort -k4,4 -k8,8nr | uniq > ${outputdir}/${outFn}
rename_uniqueRows ${outputdir}/${outFn}
