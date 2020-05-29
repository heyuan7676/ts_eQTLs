#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00

source ./GLOBAL_VAR.sh

feature=bedfiles_NullMotifs_1X_GC_Repeats_unmatched

for chr in {1..22}
do
echo chr${chr}
outputdir=${datadir}/annotations/allPairs/tissue_${feature}


#### format the final matrix

cd ${outputdir}
rm -f temp_*chr${chr}_${feature}


# header line
outFn=SNP_loc_chr${chr}_${feature}.bed
echo ${outFn}

command=paste
echo "SNP" > temp_name_chr${chr}_${feature}
for x in `ls *chr${chr}_*unique`
do
	N=`wc -l ${x} | awk '{print $1}'`
	if [[ "$N" -eq 0 ]]
	then
		echo $x
		continue	
	fi
	command="$command <(cat $x | cut -f8)"

        y=${x#SNP_loc_}
        y=${y%.bed_unique}
        echo $y >> temp_name_chr${chr}_${feature}
done

paste -sd"	" temp_name_chr${chr}_${feature} > ${outFn}

eval $command > temp_batch_chr${chr}_${feature}

paste <(awk '{print $4}' ${x}) <(cat temp_batch_chr${chr}_${feature}) >> ${outFn}
#cat ${outFn} | sed 's/\. /0      /g' > ${outFn}_temp
awk '{print NF}' ${outFn} | uniq




done
