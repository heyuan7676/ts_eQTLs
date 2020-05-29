#!/bin/bash
#SBATCH --time 30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue

### aggregate pairs for credible sets
source ./GLOBAL_VAR.sh

tissue="$1"
maindir=${datadir}/caviar_output_GTEx_LD
inputdir=${maindir}/${tissue}.allpairs.txt

function extract_pairs {
	chr="$1"
	cd ${inputdir}/chr${chr}/set
        for geneSetFn in `ls *`
        do
                geneName=${geneSetFn/_new.out_set/}
                while read -r line
                do
                        SNP=chr${chr}_"$line"
                        echo ${geneName}"	"${SNP} >> ${tissueOutFn}
                done < ${geneSetFn}
        done
}

function SNPMAF {
        tissue="$1"
        cisFn=${cisDir}/${tissue}.allpairs.txt
        paste <(cat ${cisFn} | cut -d'_' -f 1,2) <(cat ${cisFn} | cut -d'	' -f3-) | awk '{if($6>0.01) print $0}' > ${cisFn_matched_SNPs}/${tissue}.allpairs.maf0.01.txt
        cisFn=${cisFn_matched_SNPs}/${tissue}.allpairs.maf0.01txt

        outFn=${maindir}/snpMAF_filtered/${tissue}.tissuecbset.SNPMAF0.01.txt
        awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' ${tissueOutFn} ${cisFn} | sort -u  | sed 's/ /	/g' > ${outFn}
}

outputdir=${maindir}/aggregate
tissueOutFn=${outputdir}/${tissue}_95set_pairs.txt

mkdir -p ${outputdir}
rm -f ${tissueOutFn}
for chr in {1..22}
do
	extract_pairs ${chr}
done


wc -l ${tissueOutFn}

cisDir=/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations
cisFn_matched_SNPs=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/allPairs

#SNPMAF ${tissue}


