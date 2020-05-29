#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -p quickdebug

module load bedtools

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/ENCODE
tooldir=/home-4/yhe23@jhu.edu/work/yuan/tools/liftOver
outdir=${datadir}/formatting

rm -r -f ${outdir}
mkdir ${outdir}
### lift GRCh38 to hg19

cd ${tooldir}
filename=`ls ${datadir}/download | grep "bed$" `

for f in $filename
do
	echo $f
        oldfn=${datadir}/download/$f
        #cp $oldfn ${datadir}/save
        #gunzip ${oldfn}

	#oldfn=${oldfn/.gz/}
	#sed -i 's/\.0//g' ${oldfn}
	awk '{print $1,$2,$3,$4,$5}' ${oldfn} | sed 's/ /	/g' > ${oldfn}_correct

        expr=`basename ${oldfn%.bed}`
        tissue=`cat ${datadir}/metadata.tsv | grep ^${expr} | cut -d'	' -f 7`
	histone=`cat ${datadir}/metadata.tsv | grep ^${expr} | cut -d'	' -f 19`
        genome=`cat ${datadir}/metadata.tsv | grep ^${expr} | cut -d'	' -f 44`
        newfn=${outdir}/${tissue}_${genome}_${histone}_${expr}.bed

        if [[ $genome == "hg19" ]]
        then
                echo $newfn
                ./liftOver ${oldfn}_correct hg19ToHg38.over.chain.gz ${newfn//hg19/GRCh38} ${newfn%.bed}_unlifted.bed
	else
		mv ${oldfn}_correct ${newfn}
        fi
done


cd ${datadir}/formatting
mkdir -p hg19_notused
mkdir -p GRCh38


mv *GRCh38*bed GRCh38
mv *hg19*bed hg19_notused



### merge data for the same tissue

cd ${datadir}/formatting/GRCh38

rm -f prefix.txt

for x in *bed
do
        echo ${x} | awk -F '-human' '{print $1}' >> prefix.txt
done

cat prefix.txt | sort | uniq > temp
mv temp prefix.txt

while IFS= read -r line
do
        t=$line
        echo $t
        cat ${t}*bed | sort -k1,1 -k2,2n | bedtools merge -c 5 -o mean > ${t/_GRCh38/}_merged.bed
done < prefix.txt


rm -r -f ${datadir}/merged
mkdir -p ${datadir}/merged
cp ${datadir}/formatting/GRCh38/*merged.bed ${datadir}/merged


