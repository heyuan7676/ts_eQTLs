#!/bin/bash

source ./GLOBAL_VAR.sh

for g in {0..22}
do
	sbatch 6_gsea_gene.sh group${g}
done
