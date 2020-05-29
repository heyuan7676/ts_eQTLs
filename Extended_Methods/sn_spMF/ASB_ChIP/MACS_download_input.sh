#!/bin/bash
ml samtools
cd /work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/TFBS_ChIP_seq/MACS_peaks/input/

#### HNF4A

### 37 adult
wget https://www.encodeproject.org/files/ENCFF479OYS/@@download/ENCFF479OYS.bam
wget https://www.encodeproject.org/files/ENCFF864EDI/@@download/ENCFF864EDI.bam
# control
wget https://www.encodeproject.org/files/ENCFF406XEH/@@download/ENCFF406XEH.bam
mv ENCFF406XEH.bam HNF4A_S1_control.bam

samtools merge HNF4A_S1.bam ENCFF479OYS.bam ENCFF864EDI.bam

## 4 child
wget https://www.encodeproject.org/files/ENCFF658AIS/@@download/ENCFF658AIS.bam
wget https://www.encodeproject.org/files/ENCFF058OZQ/@@download/ENCFF058OZQ.bam
# control
wget https://www.encodeproject.org/files/ENCFF337WHA/@@download/ENCFF337WHA.bam
mv ENCFF337WHA.bam HNF4A_S2_control.bam

samtools merge HNF4A_S2.bam ENCFF658AIS.bam ENCFF058OZQ.bam


#### CTCF
wget https://www.encodeproject.org/files/ENCFF196CKN/@@download/ENCFF196CKN.bam
mv ENCFF196CKN.bam CTCF_S1.bam
# control
wget https://www.encodeproject.org/files/ENCFF621CXJ/@@download/ENCFF621CXJ.bam
mv ENCFF621CXJ.bam CTCF_S1_control.bam



