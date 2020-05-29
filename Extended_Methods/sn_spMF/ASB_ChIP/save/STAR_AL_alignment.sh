#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=50:00:00



## step1
#cd /home-4/yhe23@jhu.edu/work/yuan/tools/STAR-2.7.1a/bin/Linux_x86_64_static
#genome_fa_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/STAR/Homo_sapiens.GRCh38.dna.toplevel.fa
#gtf_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/STAR/gencode_chromosomenames.v30.annotation.gtf
#genome_index_idr=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/STAR/genome_index

#./STAR  --runMode genomeGenerate --runThreadN 6  --genomeDir ${genome_index_idr} --genomeFastaFiles ${genome_fa_dir} --sjdbGTFfile ${gtf_dir} --limitGenomeGenerateRAM 115999096192 

## step2
genome_index_idr=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/STAR/genome_index
vcfFn=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/STAR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_nochromsome.vcf
fastq_dir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/fastq/

cd ${fastq_dir}
fn="$1"  # fn=CTCF_Thyroid.fastq.files.txt
fastq_files=`cat ${fn}`

outputdir=/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/STAR_output
outputfn=${outputdir}/${fn%.fastq.file.txt}_ChIP_seq

cd /home-4/yhe23@jhu.edu/work/yuan/tools/STAR-2.7.1a/bin/Linux_x86_64_static
./STAR --genomeDir ${genome_index_idr}  --runThreadN 24 --readFilesIn ${fastq_files} --readFilesCommand zcat --outFileNamePrefix ${outputfn} --varVCFfile ${vcfFn} --outSAMattributes NH HI AS nM NM MD vA vG vW --outSAMtype BAM SortedByCoordinate --waspOutputMode SAMtag










