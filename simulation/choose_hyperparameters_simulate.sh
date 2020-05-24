#!/bin/bash
#SBATCH --time 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p skylake


tau="$1"
seed="$2"


xfn=simulation/input/Input_tau${tau}_seed${seed}_X.txt
wfn=simulation/input/Input_tau${tau}_seed${seed}_W.txt
factor_dir=simulation/output/tau${tau}_seed${seed}

mkdir -p ${factor_dir}

for K in {6..9}
do
for a1 in {1..10}
do
	a1=`echo "scale=1; $a1 / 10"  | bc`
        for l1 in {1..10}
        do
		l1=`echo "scale=1; $l1 / 10"  | bc`
		echo $xfn, $a1,$l1
                Rscript sn_spMF/run_MF.R -x ${xfn} -w ${wfn} -O ${factor_dir} -k ${K} -a ${a1} -l ${l1} -p 0 -c 1
        done
done
done


