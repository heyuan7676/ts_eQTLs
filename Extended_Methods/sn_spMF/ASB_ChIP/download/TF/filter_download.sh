#!/bin/bash

metafn=metadata.tsv
output=download_filtered.txt

## Hand curation
#bash change_tissue_names.sh ${metafn}


rm -f ${output}
rm -f ${output_temp}

{
read
while IFS= read line
do
        fn=`basename "$line" | cut -d. -f1`
	echo $fn

	id=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 1`
        filetype=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 2`

	tissue=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 7`
	paired=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 35`

	tf=`cat ${metafn} | grep ^${fn} | cut -d'	' -f 19`

	if [[ "$filetype" == "fastq" ]]
	then
                echo $line >> ${output}
		echo ${id}, ${tf}, ${tissue}, ${filetype}, ${paired}
        fi

done
} < download.txt


