#!/bin/bash
#$ -N cuffdiff
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load cufflinks

cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf C1_R1_thout/accepted_hits.bam, C2_R1_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf C1_R1_thout/accepted_hits.bam, C2_R1_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf C1_R1_thout/accepted_hits.bam, C2_R1_thout/accepted_hits.bam
