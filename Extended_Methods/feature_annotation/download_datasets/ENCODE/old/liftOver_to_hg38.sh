#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=50:00:00

feature="$1"
featuredir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v7/datasets/annotations/${feature}
outdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/${feature}
tooldir=/home-4/yhe23@jhu.edu/work/yuan/tools/liftOver

rm -f -r ${outdir}
mkdir -p ${outdir}

cd ${featuredir}
fns=`ls *bed`

for f in ${fns}
do
	echo $f
	cd ${tooldir}
	./liftOver ${featuredir}/${f} hg19ToHg38.over.chain.gz ${outdir}/${f} ${outdir}/${f%.bed}_unlifted.bed
	cd ${outdir}
	sort -k1,1 -k2,2n ${f} > tmp_${f}
	mv tmp_${f} ${f}
done


cd ${outdir}
mkdir -p unlifted
mv *unlifted* unlifted/
