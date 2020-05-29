#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

ml parallel

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif
rm -f formatTF_batch_commands.txt

#for f in `ls ${datadir}/tsv/*tsv.gz`
for f in `cat makeup_TFs.txt`
do
	f=`echo ${f} | cut -d'_' -f1`
	f=${f}.tsv.gz
	f=`basename ${f}`
        TFname=`zcat ${datadir}/tsv/$f | head -n2 | tail -n1 | awk '{print $4}'`
        newf=${f%.tsv.gz}_${TFname}.bed
        newf=${newf/(var.2)/-var.2}
        newf=${newf/(var.3)/-var.3}
	echo "zcat ${datadir}/tsv/$f | sed "1d" | cut -f1,2,3,6 | sort -k1,1 -k2,2n > ${datadir}/bed/${newf}" >> formatTF_batch_commands.txt
done


parallel :::: formatTF_batch_commands.txt
