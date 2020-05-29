#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH -p lrgmem

source ./GLOBAL_VAR.sh

feature="$1"
subfeature="$2"
inputfn="$3"
dir="$4"
comman_string=${inputfn%.bed}

featuredir=${datadir}/annotations/${feature}/merged
inputdir=${datadir}/prediction_model/${dir}

outputdir=${datadir}/annotations/prediction_model/tissue_${feature}
rm -f ${outputdir}/${comman_string}*${subfeature}*

mkdir -p ${outputdir}

cd ${featuredir}
for featurefn in `ls *${subfeature}*bed`
do
        #echo $featurefn
	outFn=${comman_string}_${featurefn}
	/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wao -a ${inputdir}/${inputfn} -b ${featurefn} -nonamecheck -sorted  | sed 's/\.	/0	/g' > ${outputdir}/${outFn}
	done



#### format the matrix of peaks 
cd ${outputdir}


command=paste
namefn=temp_name_${comman_string}_${subfeature}
echo "SNP" > ${namefn}

for fn in `ls ${comman_string}*${subfeature}_merged.bed`
do
        #echo $fn
	#awk '!seen[$4]++' $fn | awk '{print $8";"$4}'> temp_$fn
        awk '$4 != prev {if (NR != 1) print ";"prev; prev=$4; delete a};
                !($8 in a){a[$8]++; printf "%s ", $8};
                END {print ";"prev}' ${fn} | sed 's/ /,/g'> tp_${fn}
	# content
	command="$command <(cat tp_$fn | cut -d';' -f1)"
	# colnames
	y=${fn#${comman_string}_}
        y=${y%.bed}
        y=${y%_merged}
        echo $y >> ${namefn}

done

# header line
outFn=${comman_string}_${subfeature}.bed
paste -sd"      " ${namefn} > ${outFn}

bodyfn=temp_batch_${comman_string}_${subfeature}
eval $command > ${bodyfn}
paste <(cat tp_${fn} | cut -d';' -f2) <(cat ${bodyfn}) >> ${outFn}

cat ${outFn} | sed 's/ /	/g' > tp_${comman_string}
mv tp_${comman_string} ${outFn}

awk '{print NF}' ${outFn} | uniq

rm -f ${namefn}
rm -f ${bodyfn}
rm -f tp_${comman_string}*${subfeature}_merged.bed


### format the matrix of coverage
command=paste
namefn=temp_name_${comman_string}_${subfeature}
echo "SNP" > ${namefn}

for fn in `ls ${comman_string}*${subfeature}_merged.bed`
do
        #echo $fn
        #awk '!seen[$4]++' $fn | awk '{print $8";"$4}'> temp_$fn
        awk '$4 != prev {if (NR != 1) print ";"prev; prev=$4; delete a};
                !($9 in a){a[$9]++; printf "%s ", $9};
                END {print ";"prev}' ${fn} | sed 's/ /,/g'> tp_${fn}
        # content
        command="$command <(cat tp_$fn | cut -d';' -f1)"
        # colnames
        y=${fn#${comman_string}_}
        y=${y%.bed}
        y=${y%_merged}
        echo $y >> ${namefn}

done

# header line
outFn=${comman_string}_${subfeature}_coverage.bed
paste -sd"	" ${namefn} > ${outFn}

bodyfn=temp_batch_${comman_string}_${subfeature}
eval $command > ${bodyfn}
paste <(cat tp_${fn} | cut -d';' -f2) <(cat ${bodyfn}) >> ${outFn}

cat ${outFn} | sed 's/ /	/g' > tp_${comman_string}
mv tp_${comman_string} ${outFn}
awk '{print NF}' ${outFn} | uniq

rm -f ${namefn}
rm -f ${bodyfn}
rm -f tp_${comman_string}*${subfeature}_merged.bed
mv ${comman_string}_*_${subfeature}_merged.bed each_tissue


