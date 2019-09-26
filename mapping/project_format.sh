#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH -p lrgmem

cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL/

FM_fn=SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110.mat
head -n1 ${FM_fn%.mat}_startidx20_LD1_Loadings_projection.txt > ${FM_fn%.mat}_LD1_Loadings_projection.txt
for x in ${FM_fn/.mat/}_startidx*_LD1_Loadings_projection.txt
do
	sed 1d ${x} >> ${FM_fn%.mat}_LD1_Loadings_projection.txt
done

cat ${FM_fn%.mat}_LD1_Loadings_projection.txt | sed 's/ /	/g' | sed 's/"//g' > temp
mv temp ${FM_fn%.mat}_LD1_Loadings_projection.txt 

mv ${FM_fn/.mat/}_startidx*_LD1_Loadings_projection.txt startidx/
