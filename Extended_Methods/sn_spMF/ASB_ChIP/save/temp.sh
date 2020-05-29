#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

alignmetn_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output


ml python/2.7
for fn in `ls ${alignmetn_dir}/*.filtered.txt | xargs -n1 basename`
do
	echo $fn
        python 5_TFBS_ChIP_ASB_format_samFiles.py ${fn}
done

