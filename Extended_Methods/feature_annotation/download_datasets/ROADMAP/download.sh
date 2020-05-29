#!/bin/bash
#SBATCH -p lrgmem


celltype="$1"
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/${celltype}_18_core_K27ac_hg38lift_mnemonics.bed.gz

