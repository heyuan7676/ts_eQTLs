#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1


source ./GLOBAL_VAR.sh
#LMfn=flashr_Loadings_beta_BH_alpha0.05_corrected
#LMfn=Thresholding_eGenes
group="$1"

geneset=topGene5
#geneset="$2"
#LMfn="$3"
#tis="$1"
gsea_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/
gsea_fn=c5.bp.v6.2.symbols.gmt.txt

file1=${LMfn}_${geneset}_${group}.txt
#file1=Thresholding_eGenes_topGene30_${tis}.txt
file3=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/tested_genes_in_cbset/Group${group/group/}_genes.txt

outdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/GSEA
outfn=${outdir}/${file1%.txt}_enrichment_${gsea_fn}

echo $outfn

module load python/2.7
python run_enrichment.py ${pairdir}/${file1} ${file3} ${gsea_dir}/${gsea_fn} ${outfn}

cat ${outfn} | head -n5 | cut -d'	' -f1-7,9
