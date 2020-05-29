#!/bin/bash


metafn=metadata.tsv
output=metadata_filtered.txt
output_temp=metadata_filtered_out.txt
download_fn=download_filtered.txt


## Hand curation
#bash change_tissue_names.sh ${metafn}


rm -f ${output}
rm -f ${output_temp}
rm -f ${download_fn}


cat ${metafn} | grep "bed narrowPeak" > temp
mv temp ${metafn}

{
read
while IFS= read -r line
do
        fn=`basename "$line" | cut -d. -f1`
	echo $fn

	id=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 1`
        filetype=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 3`
	filestatus=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 48`
	info=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 52`
	errorInfo=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 53`
	tissue=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 7`
	histone=`cat ${metafn} | grep ${fn} | cut -d'	' -f 19`
	downloadaddress=`cat ${metafn} | grep ${fn} | cut -d'	' -f 43`

	if [[ "$filetype" == "stable peaks" ]]  && [[ "$filestatus" == "released" ]] && [[ ${errorInfo} != *"extremely"* ]]
        then
                echo $line >> ${output}
		echo ${downloadaddress} >> ${download_fn}
	else
		echo ${id}, $tissue, $filetype, $filestatus, $errorInfo >> ${output_temp}
        fi

done
} < ${metafn}

