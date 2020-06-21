#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -p shared


FM_fn="$1"
echo ${FM_fn}
source ./GLOBAL_VAR.sh

### fit the loadings
cd ${scripts_dir}/../${FM_fn} 
command_fn=fitL_commands.txt
rm -f ${command_fn}
for i in {0..36}
do
	outfn=${ll_dir}/${FM_fn}_startidx${i}_Loadings_beta.txt
	if [ ! -f ${outfn} ]
	then
		echo "Rscript ${scripts_dir}/0_loadings/lm.R ${i} ${FM_fn}" >> ${command_fn}
	fi
done

ml R/3.5.1
ml parallel
parallel -j 3 :::: ${command_fn}


### format the loadings
cd ${ll_dir}
for feature in pvalue beta
do
        head -n1 ${FM_fn}_startidx0_Loadings_${feature}.txt > ${FM_fn}_Loadings_${feature}.txt
	for i in {0..36}
        do
		x=${FM_fn}_startidx${i}_Loadings_${feature}.txt
                sed 1d ${x} >> ${FM_fn}_Loadings_${feature}.txt
        done
        mv  ${FM_fn}_startidx*_Loadings_${feature}.txt startidx/
done


### BH correction
Rscript ${scripts_dir}/0_loadings/adj_p.R 0.05 ${FM_fn}









