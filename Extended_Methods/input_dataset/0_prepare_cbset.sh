#!/bin/bash
#SBATCH --time 30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --no-requeue

source ./GLOBAL_VAR.sh
NBlocks=38
ml python/2.7
ml parallel


### Step 1)
### aggregate pairs for credible sets

fn=step1_commands.txt
rm -r ${fn}
for t in `cat tissues.txt`
do
	echo "bash step1_cbset_each_tissue.sh ${t}"  >> ${fn}
done
parallel :::: ${fn}


### Step 2)
### Filter genes to keep only protein coding genes
bash step2.1_aggreate_filterGenes_noMAF.sh


### Step 2.2)
### collect cis-eQTL results for all eQTL pairs in the credible set
### parallel implementation in the script
bash step2_collect_cisresults.sh


### Step 3)
### Merge info across tissues
bash step3_merge_files.sh


### Step 4) -- deal with SNPs in LD blocks
### Step 4.1)
### Compute R2 for SNPs
r2_thr=0.2
fn=SNPLD_commands.txt
rm -r ${fn}
for chr in {1..22}
do
	echo "bash SNP_LD.sh ${r2_thr} ${chr}" >> ${fn}
done
parallel :::: ${fn}


genotypedir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes
for r2 in 1
do
	# collect the results
        summary_fn=${genotypedir}/plink_r2_${r2Thr}.ld
        rm -r ${summary_fn}
        for chr in {1..22}
        do
                r2Out=${genotypedir}/plink_r2_${r2Thr}_bychr/plink_r2_${r2Thr}_${chr}.ld
                sed 1d ${r2Out} | awk -v var="$r2" '{if($7>var) print $3,$6,$7}'  >> ${summary_fn}
        done

	### Step 4.2)
	### Compute connected components for SNPs
	python SNP_LD_comp.py ${r2}

	### Step 4.3)
	### For each LD block, use the pair with minimum genometric mean p-value, and save the others
	outdir=${datadir}/input_pairs_fitModel
	mkdir -p ${outdir}

	fn=filter_SNP_groupSNPs_${r2}_commands.txt
	rm -f ${fn}
	for ith in {0..NBlocks}
	do
		echo "python filter_SNP_groupSNPs.py ${r2} ${suffix} ${ith}" >> ${fn}
	done
	parallel :::: ${fn}

	### Step 4.4)
	### Construct input data for factor learning, and for loading learning
	python filter_SNPs_TagSNPs.py ${r2} ${suffix}
done



### Step 5
### Extract AF information for tested SNPs
bash SNP_MAF.sh
