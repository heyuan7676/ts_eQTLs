#!/bin/bash


metafn=metadata.tsv
metafn_new=${metafn%.tsv}_final.tsv
output_temp=metadata_filtered_out.txt
output=file.txt

## Hand curation
#bash change_tissue_names.sh ${metafn}


rm -f ${output}
rm -f ${output_temp}
rm -f ${metafn_new}
rm -f ${output}

{
read
while IFS= read -r line
do
        fn=`basename "$line" | cut -d. -f1`
	echo $fn

	id=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 1`
        filetype=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 3`
	info=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 52`
	errorInfo=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 53`
	tissue=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 7`
	histone=`cat ${metafn} | grep ${fn} | cut -d'	' -f 19`
	downloadaddress=`cat ${metafn} | grep ${fn} | cut -d'	' -f 43`

	if [[ "$filetype" == "peaks" ]] && [[ ${errorInfo} != *"extremely"* ]]
        then
		echo ${id}, $tissue, $filetype, $errorInfo >> ${metafn_new}
		echo ${downloadaddress} >> ${output}
	else
		echo ${id}, $tissue, $filetype, $errorInfo >> ${output_temp}
        fi

done
} < ${metafn}

