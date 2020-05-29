#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -p lrgmem

### produce:  SNPs x Tissues in DNase data


ml parallel
source ./GLOBAL_VAR.sh

command_fn=extractSNP_features_commands.txt
rm -f ${command_fn}

anno_feature="$1"
datadir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/
annodir=${datadir}/annotations/allPairs/TFBS_ChIP
annoFn=SNP_loc_${anno_feature}.bed


iddir=${pairdir}
random="$2"
echo "extract SNP features"
for group in {0..22}
do
		echo "group"$group
		idfn=${LMfn}_outlierPairs${random}_group${group}.txt
		outfn=${LMfn}_outlierSNPs${random}_${anno_feature}_group${group}.txt
		rm -f ${outfn}
		bash 4_extractSNP_feature.sh ${iddir} ${idfn} ${annodir} ${annoFn} ${sigSNPfeaturedir} ${outfn}

done

