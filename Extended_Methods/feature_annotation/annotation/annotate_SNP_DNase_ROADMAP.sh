#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00

##### ROADMAP annotataion

#http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
#http://egg2.wustl.edu/roadmap/figures/mainFigs/Figure_2.jpg

#for x in *gz; do gunzip $x; done

source ./GLOBAL_VAR.sh

feature=bedfiles_NullEnh_1X
for chr in {1..22}
do

featuredir=${datadir}/annotations/${feature}

inputdir=${datadir}/GTEx_datasets/allPairs/maf_dis/SNP_split_by_chr/
inputfn=SNP_loc_chr${chr}.bed

outputdir=${datadir}/annotations/allPairs/tissue_${feature}
mkdir -p ${outputdir}
rm -f ${outputdir}/*chr${chr}_*
mkdir -p ${outputdir}

cd ${featuredir}
for featurefn in `ls *bed`
do
        #echo $featurefn
	#sort -k1,1 -k2,2n $featurefn > tmp_$featurefn
	#mv tmp_$featurefn $featurefn
	outFn=${inputfn%.bed}_${featurefn}
	/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wao -a ${inputdir}/${inputfn} -b ${featurefn} -nonamecheck -sorted  | sed 's/\.	/0	/g' > ${outputdir}/${outFn}
done



#### format the final matrix
cd ${outputdir}
rm -f temp_*chr${chr}_*

for fn in `ls *chr${chr}_*bed`
do
        echo $fn
	#awk '!seen[$4]++' $fn | awk '{print $8";"$4}'> temp_$fn
        awk '$4 != prev {if (NR != 1) print ";"prev; prev=$4; delete a};
                !($8 in a){a[$8]++; printf "%s ", $8};
                END {print ";"prev}' ${fn} | sed 's/ /,/g'> temp_${fn}
done

command=paste
for i in temp*chr${chr}_*bed
do 
       command="$command <(cat $i | cut -d';' -f1)"
done
eval $command > temp_batch_chr${chr}_${feature}


# header line
outFn=${inputfn/.bed/}_${feature}.bed

echo "SNP" > temp_name_chr${chr}_${feature}
for x in temp_SNP_loc_chr${chr}_*.bed
do
        #echo $x
        y=${x#temp_SNP_loc_}
        y=${y%.bed}
        y=${y%_merged}
        echo $y >> temp_name_chr${chr}_${feature}
done
paste -sd"      " temp_name_chr${chr}_${feature} > ${outFn}

paste <(cat temp_${fn} | cut -d';' -f2) <(cat temp_batch_chr${chr}_${feature}) >> ${outFn}
cat ${outFn} | sed 's/ /	/g' > temp_chr${chr}
mv temp_chr${chr} ${outFn}
awk '{print NF}' ${outFn} | uniq

rm -f temp_*chr${chr}_*
done
