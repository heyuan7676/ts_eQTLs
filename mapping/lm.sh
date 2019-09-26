#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3

ml gcc/5.5.0
ml R/3.5.1


#FM_fn=SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_topPair_K25_a125_l15000.mat
FM_fn =SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110.mat
inputprefix=v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks
r2=1

idx="$1"
Rscript lm.R ${FM_fn} ${inputprefix} ${r2} ${idx}









