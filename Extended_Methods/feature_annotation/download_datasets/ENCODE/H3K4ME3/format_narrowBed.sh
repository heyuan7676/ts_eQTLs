#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -p quickdebug

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/ENCODE
metafn=H3K4me3_metadata_final_GM12878.tsv

filename=`ls ${datadir}/download | grep "bed$" `

for f in $filename
do
	echo $f
        oldfn=${datadir}/download/$f

	awk '{print $1,$2,$3,$4,$9}' ${oldfn} | sed 's/ /	/g' > ${oldfn}_correct

        expr=`basename ${oldfn%.bed}`
        tissue=`cat ${metafn} | grep ^${expr} | cut -d'	' -f 7`
	tissue=${tissue// /}
	tissue=${tissue//\'/}
	echo $tissue
	expi=`cat ${metafn} | grep ^${expr} | cut -d'	' -f 1`
        newfn=${datadir}/H3K4ME3_${tissue}_${expi}.bed

	sort -k1,1 -k2,2n ${oldfn}_correct > ${newfn}
done

rm ${datadir}/download/*_correct
