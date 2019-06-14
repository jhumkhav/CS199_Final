#!/bin/bash
#$ -N cufflinks
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load cufflinks

cufflinks -p 8 -o C1_R1_clout C1_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R2_clout C1_R2_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R3_clout C1_R3_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R1_clout C2_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R2_clout C2_R2_thout/accepted_hits.bam

