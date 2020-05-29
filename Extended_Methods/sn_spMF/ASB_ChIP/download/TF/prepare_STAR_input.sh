#!/bin/bash

metafn=metadata.tsv
suffix_fn=.fastq.files.txt

dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/fastq/

rm -f *${suffix_fn}


{
read
while IFS= read line
do
        fn=`basename "$line" | cut -d. -f1`
	echo $fn

	tissue=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 7`
	paired=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 35`
	tf=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 19 | cut -d'-' -f 1`

	if [[ "$paired" == "paired-ended" ]]
	then
		other_fn=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 37`
		echo ${dir}${fn}.fastq.gz" "${dir}${other_fn}.fastq.gz >> ${tf}_${tissue}${suffix_fn}
	else
		echo ${dir}${fn}.fastq.gz >> ${tf}_${tissue}${suffix_fn}
        fi

done
} < download_filtered.txt

for f in *${suffix_fn}
do
	temp=`cat ${f} | tr '\n' ',' `
	line=`echo $temp | rev | cut -c 2- | rev`
	echo ${line} > ${f}
done
