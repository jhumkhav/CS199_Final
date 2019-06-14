#!/bin/bash
#$ -N samtools
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load samtools

samtools sort -@ 8 -o C1_R1.bam C1_R1.sam
samtools sort -@ 8 -o C1_R2.bam C1_R2.sam
samtools sort -@ 8 -o C1_R3.bam C1_R3.sam
samtools sort -@ 8 -o C2_R1.bam C2_R1.sam
samtools sort -@ 8 -o C2_R2.bam C2_R2.sam
samtools sort -@ 8 -o C2_R3.bam C2_R3.sam
