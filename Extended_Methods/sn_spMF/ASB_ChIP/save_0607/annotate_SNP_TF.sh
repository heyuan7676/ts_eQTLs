#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/
featuredir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq

TFfn="$1"
for chr in {1..22}
do
	echo chr${chr}

	inputdir=${datadir}/GTEx_datasets/allPairs/maf_dis/SNP_split_by_chr/
	inputfn=SNP_loc_chr${chr}.bed

	outputdir=${datadir}/annotations/allPairs/TFBS_ChIP
	outFn=${inputfn%.bed}_${TFfn}

	/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wo -a ${inputdir}/${inputfn} -b ${featuredir}/${TFfn} -nonamecheck -sorted | sort -k4,4 -k8,8nr | uniq > ${outputdir}/${outFn}

done


cd ${outputdir}
cat SNP_loc_chr*_${TFfn} > SNP_loc_${TFfn}
mv SNP_loc_chr*_${TFfn} split_by_chr/
