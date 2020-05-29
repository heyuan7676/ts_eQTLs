#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1


source ./GLOBAL_VAR.sh
LMfn=Thresholding_eGenes
tis="$1"
gsea_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/
gsea_fn=c5.bp.v6.2.symbols.gmt.txt

file1=ts_closeToTop_FOLDS100_PROP0.5_PVALUE0.001_N15_topGene4_group${tis}.txt
file3=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs/tested_genes_in_Thresholding/${tis}_genes.txt

outdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/enrichmentTest/GSEA
outfn=${outdir}/${file1%.txt}_enrichment_${gsea_fn}

echo $outfn

module load python/2.7
python run_enrichment.py ${pairdir}/${file1} ${file3} ${gsea_dir}/${gsea_fn} ${outfn}

cat ${outfn} | head -n5 | cut -d'	' -f1-7,9



