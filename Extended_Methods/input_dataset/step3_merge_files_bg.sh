#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --no-requeue
#SBATCH -p debug

source ./GLOBAL_VAR.sh

inputdir=${datadir}/cbset_ciseQTL_results
outdir=${datadir}/input_pairs
rm -r -f ${outdir}
mkdir -p ${outdir}


### distance to TSS
awk '{print $1, $2, $3}' ${inputdir}/Adipose_Subcutaneous.${suffix}.txt | sed 's/ /	/g' >  ${outdir}/${suffix}.tss_distance.txt


### get the pairs

pairsFn=${inputdir}/Adipose_Subcutaneous.${suffix}.txt
command=paste
command="$command <(sed "1d" $pairsFn | cut -f1,2)"


### merge cols from multiple files

command_pvalue=${command}
command_slope=${command}
command_se=${command}

tissues=`cat tissues.txt`
for t in ${tissues}
do
        i=${inputdir}/${t}.${suffix}.txt
        command_pvalue="$command_pvalue <(sed "1d" $i | cut -f7)"
	command_slope="$command_slope <(sed "1d" $i | cut -f8)"
	command_se="$command_se <(sed "1d" $i | cut -f9)"
done

eval $command_pvalue | sed 's/ /	/g' > ${outdir}/${suffix}.pvalue.txt
eval $command_slope  | sed 's/ /	/g' > ${outdir}/${suffix}.slope.txt
eval $command_se     | sed 's/ /	/g' > ${outdir}/${suffix}.se.txt


