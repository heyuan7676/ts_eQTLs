#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH -p debug




source ./GLOBAL_VAR.sh
source ./utils.sh

idsDir=${datadir}/caviar_output_GTEx_LD/aggregate
idsFn=${idsDir}/v8_cbset_95_allPairs_filteredGenes.txt

cisDir=/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations
cisFn_matched_SNPs=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/allPairs

outDir=${datadir}/cbset_ciseQTL_results
mkdir -p ${outDir}

ml parallel
rm -f step2_collect_cisresults_commands.txt
for t in `cat tissues.txt`
do
	echo "indexPairsFromciseQTL ${t}" >> step2_collect_cisresults_commands.txt
done
parallel :::: step2_collect_cisresults_commands.txt



