#!/bin/bash
#SBATCH --time 50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

### produce:  SNPs x Tissues in DNase data


ASB_SNP_fn="$1"
specifideQTL_inGWAS=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/gwasSNP
for fn in `ls ${specifideQTL_inGWAS}/group15* | grep -v random | xargs -n1 basename`
do
	grep -Ff ${ASB_SNP_fn} ${specifideQTL_inGWAS}/${fn} > intersection_result/${ASB_SNP_fn}_${fn}
	N=`wc -l intersection_result/${ASB_SNP_fn}_${fn} | awk '{print $1}'`
	if [ "$N" -eq "0" ]
	then
		rm intersection_result/${ASB_SNP_fn}_${fn}
	fi
done

