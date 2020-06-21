#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p shared


FM_fn="$1"
i="$1"
echo ${FM_fn}
source ./GLOBAL_VAR.sh

### fit the loadings
cd ${scripts_dir}/../${FM_fn} 

ml R/3.5.1
outfn=${ll_dir}/${FM_fn}_startidx${i}_Loadings_beta.txt
if [ ! -f ${outfn} ]
then
	Rscript ${scripts_dir}/0_loadings/lm.R ${i} ${FM_fn}
fi

