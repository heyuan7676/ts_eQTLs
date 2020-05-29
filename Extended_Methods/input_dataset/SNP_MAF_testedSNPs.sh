#!/bin/bash

iddir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs
idfn=v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete.pvalue.geman.txt
inputdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes/vcf_chr
inputfn=v8_SNP_AF.txt 

outdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs
outfn=v8_cbset_95_SNPs_AF.txt


awk -F '\t' 'NR==FNR {id[$2]; next} $1 in id' ${iddir}/${idfn} ${inputdir}/${inputfn} > ${outdir}/${outfn}
