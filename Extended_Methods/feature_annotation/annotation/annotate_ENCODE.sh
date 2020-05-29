#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH -p lrgmem

##### ROADMAP annotataion

#http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
#http://egg2.wustl.edu/roadmap/figures/mainFigs/Figure_2.jpg

#for x in *gz; do gunzip $x; done

source ./GLOBAL_VAR.sh

histone="$1"
sett="$2"
chr="$3"

feature=ENCODE_${histone}_${sett}

featuredir=${datadir}/annotations/ENCODE/${histone}/${sett}

inputdir=${datadir}/GTEx_datasets/allPairs/maf_dis/SNP_split_by_chr/
inputfn=SNP_loc_chr${chr}.bed

outputdir=${datadir}/annotations/allPairs/tissue_${feature}
mkdir -p ${outputdir}
rm -f ${outputdir}/*chr${chr}_*
mkdir -p ${outputdir}

cd ${featuredir}
for featurefn in `ls *bed`
do
        echo $featurefn
	outFn=${inputfn%.bed}_${featurefn}
	awk '{print $1,$2,$3,$4}' ${featurefn} | sed 's/ /	/g' > temp_${chr}
	/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wao -a ${inputdir}/${inputfn} -b temp_${chr} -nonamecheck -sorted  | sed 's/\.	/0	/g' > ${outputdir}/${outFn}
	rm temp_${chr}
done



#### format the final matrix
cd ${outputdir}
rm -f temp_*chr${chr}_*

for fn in `ls *chr${chr}_*bed`
do
        #echo $fn
	#awk '!seen[$4]++' $fn | awk '{print $8";"$4}'> temp_$fn
	awk '!seen[$4]++' $fn > temp_${fn}
        #awk '$4 != prev {if (NR != 1) print ";"prev; prev=$4; delete a};
                #!($8 in a){a[$8]++; printf "%s ", $8};
                #END {print ";"prev}' ${fn} | sed 's/ /,/g'> temp_${fn}
	sed 's/	/ /g' temp_${fn} > ${fn}
	rm temp_${fn}
done

command=paste
for i in `ls *chr${chr}_*bed`
do 
       command="$command <(cat $i | cut -d'	' -f8)"
done
eval $command > temp_batch_chr${chr}_${feature}


# header line
outFn=${inputfn/.bed/}_${feature}.bed

echo "SNP" > temp_name_chr${chr}_${feature}
for x in SNP_loc_chr${chr}_*.bed
do
        #echo $x
        y=${x#SNP_loc_}
	x=`echo $x | cut -d'-' -f1`
	y=${x/SNP_loc_chr${chr}_/}
        echo $y >> temp_name_chr${chr}_${feature}
done

paste -sd"      " temp_name_chr${chr}_${feature} > ${outFn}
paste <(cat ${fn} | cut -d' ' -f4) <(cat temp_batch_chr${chr}_${feature}) | sed 's/ /	/g' >> ${outFn}


echo $outFn
awk '{print NF}' ${outFn} | uniq
#rm -f temp_*chr${chr}_*

