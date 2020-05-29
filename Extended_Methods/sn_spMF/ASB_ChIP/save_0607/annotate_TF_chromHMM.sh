#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=50:00:00
#SBATCH -p lrgmem

datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/
featuredir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq

TFfn="$1"
for tissue in `cat tissues.txt`
do
	echo ${tissue}

	inputdir=${datadir}/annotations/ROADMAP
	inputfn=${tissue}.bed

	outputdir=${featuredir}/chromHMM
	outFn=ROADMAP_${inputfn%.bed}_${TFfn}

	/home-4/yhe23@jhu.edu/work/yuan/tools/bedtools2/bin/bedtools intersect -wo -a ${inputdir}/${inputfn} -b ${featuredir}/${TFfn} -nonamecheck -sorted | sort -k4,4 -k8,8nr | uniq > ${outputdir}/${outFn}

done

