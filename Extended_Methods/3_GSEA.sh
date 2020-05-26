#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1


source ./GLOBAL_VAR.sh
module load python/2.7

tis="$1"
bg_file=${inputdatadir}/Revision/tested_genes/${tis}_genes_${FMfn}.txt


topGene_N=topGene"$2"
file1=${LMfn}_${topGene_N}_group${tis}.txt
outfn=${gsea_dir}/${file1%.txt}_enrichment.txt

echo $outfn
python ${scripts_dir}/run_enrichment.py ${pairdir}/${file1} ${bg_file} ${gsea_fn} ${outfn}
cat ${outfn} | head -n5 | cut -d'	' -f1-7,9




topGene_N=topGene30
file1=${LMfn}_${topGene_N}_group${tis}.txt
outfn=${gsea_dir}/${file1%.txt}_enrichment.txt

echo $outfn
python ${scripts_dir}/run_enrichment.py ${pairdir}/${file1} ${bg_file} ${gsea_fn} ${outfn}
cat ${outfn} | head -n5 | cut -d'	' -f1-7,9

