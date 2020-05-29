#!/bin/bash


metafn=metadata.tsv
output=metadata_filtered.txt
output_temp=metadata_filtered_out.txt
download_fn=download_filtered.txt


## Hand curation
#bash change_tissue_names.sh ${metafn}


rm -f ${output}
rm -f ${output_temp}

{
read
while IFS= read line
do
        fn=`basename "$line" | cut -d. -f1`

	id=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 1`
        filetype=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 2`
	filestatus=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 41`
	spotscore=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 42`
	errorInfo=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 43`
	tissue=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 7`
	replicates=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 25`
        sample_age=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 9`
	info=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 3`

	if [[ "$filetype" == "bed narrowPeak" ]]  && [[ "$filestatus" == "released" ]] && [[ ${errorInfo} != *"extremely"* ]]
        then
                echo $line >> ${output}
	else
		echo ${id}, $tissue, $filetype, $info, $filestatus, $errorInfo >> ${output_temp}
        fi

done
} < ${metafn}


rm -f ${download_fn}
for col in {23..29}
do
	cat metadata_filtered.txt | cut -d' ' -f ${col} | grep "download" >> ${download_fn}
done

