#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=50:00:00
#SBATCH -p shared

fn=HNF4A_Liver_single-ended_ENCFF435TOG_ENCFF366KUG_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_combined_formatted
for filter in 8 10 12
do
	bash save_reads_count.sh ${fn} ${filter}
done

fn=HNF4A_Liver_single-ended_ENCFF500ZBE_ENCFF302XOK_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_combined_formatted
for filter in 8 10 12
do
        bash save_reads_count.sh ${fn} ${filter}
done


fn=CTCF_Liver_paired-ended_ENCFF148DET_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_formatted
for filter in 8 10 12
do
        bash save_reads_count.sh ${fn} ${filter}
done


#fn=CTCF_Liver_single-ended_ENCFF002EWZ_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_formatted
#for filter in 8 10 12
#do
#        bash save_reads_count.sh ${fn} ${filter}
#done

#fn=CTCF_Liver_single-ended_ENCFF002EXB_ChIP_seqAligned.sortedByCoord.out.bam.filtered.txt_formatted
#for filter in 8 10 12
#do
#        bash save_reads_count.sh ${fn} ${filter}
#done
