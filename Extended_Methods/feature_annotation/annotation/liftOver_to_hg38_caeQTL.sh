#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00

featuredir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v7/20181116/data_new/PancanQTL
tooldir=/home-4/yhe23@jhu.edu/work/yuan/tools/liftOver

cd ${featuredir}
fns=`ls *xls`

for f in ${fns}
do
	echo $f
	awk -v s=1 '{print $2,$3,$3+s,$1}' ${f} | sed "1d" | sed 's/ /	/g' | sort -k1,1 -k2,2n > ${f/xls/bed}	
	cd ${tooldir}
	./liftOver ${featuredir}/${f/xls/bed} hg19ToHg38.over.chain.gz ${featuredir}/${f%.xls}_hg38.bed ${featuredir}/${f%.xls}_unlifted.bed
	cd ${featuredir}
done


cd ${featuredir}
mkdir -p unlifted
mv *unlifted* unlifted/
