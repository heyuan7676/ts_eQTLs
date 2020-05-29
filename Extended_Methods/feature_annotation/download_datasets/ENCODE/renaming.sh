#!/bin/bash



while IFS= read -r line
do
	t=$line
	t1=`echo $t | awk '{print $1}'`
	t2=`echo $t | awk '{print $2}'`
	cp save/H3K4ME3_${t1}_union.bed ${t2}-H3K4ME3_${t1}_union.bed
done < renaming.txt
