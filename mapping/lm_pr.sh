#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH -p lrgmem

ml gcc/5.5.0
ml R/3.5.1


FM_fn=SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110.mat
#FM_fn=SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_topPair_K25_a125_l15000.mat
inputprefix=v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks
r2=1

fn=lm_commands.txt
rm -f ${fn}
#Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 0 >> ${fn}
#for g in {1..36}
#do
echo "Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 2" >> ${fn}
echo "Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 3" >> ${fn}
echo "Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 19" >> ${fn}
echo "Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 22" >> ${fn}
echo "Rscript lm.R ${FM_fn} ${inputprefix} ${r2} 23" >> ${fn}
#done

ml parallel
parallel :::: ${fn}


cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL/
for feature in pvalue beta
do
        head -n1 ${FM_fn/.mat/}_startidx0_LD${r2}_Loadings_${feature}.txt > ${FM_fn/.mat/}_LD${r2}_Loadings_${feature}.txt
        for x in ${FM_fn/.mat/}_startidx*_LD${r2}_Loadings_${feature}.txt
        do
                sed 1d ${x} >> ${FM_fn/.mat/}_LD${r2}_Loadings_${feature}.txt
	done
	mv ${FM_fn/.mat/}_startidx*_LD${r2}_Loadings_${feature}.txt startidx/
done


cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/codes/clean_code/2_loading_learning
ml gcc/5.5.0
ml R/3.5.1

Rscript adj_p.R ${FM_fn} ${r2} 0.05
Rscript adj_p.R ${FM_fn} ${r2} 0.1
