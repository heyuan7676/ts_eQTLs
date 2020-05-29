#!/bin/bash

group="$1"

gsea_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/
gsea_fn=c5.bp.v6.2.symbols.gmt.txt

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/topTFs

#file1=tfs_gsea/group${group}_sigTFs.txt
#file1=tfs_gsea/sig_empirical_pv_liver.txt
#file1=temp.txt
#file2=tfs_gsea/group${group}_testedTFs.txt
#file2=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/tested_genes_in_cbset/Group${group/group/}_genes.txt


file1=group${group}.txt
file2=group${group}_unsig.txt

cat /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/plots/Fig5_TF_enrichment_ROADMAP_7_Enh_TPM-1_iter50.txt | awk -v var="$group" '{if($8==var) print $0}'  | awk '{if($7>=0.05) print $1,$1}' | sed 's/ /	/g'> ${file2}
cat /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/plots/Fig5_TF_enrichment_ROADMAP_7_Enh_TPM-1_iter50.txt | awk -v var="$group" '{if($8==var) print $0}'  | awk '{if($7<0.05) print $1,$1}' | sed 's/ /	/g' > ${file1}

outfn=tfs_gsea/gsea_group${group}.txt

module load python/2.7
python run_enrichment.py ${file1} ${file2} ${gsea_dir}/${gsea_fn} ${outfn}

echo ${outfn}
cat ${outfn} | head -n5 | cut -d'	' -f1-7,9

