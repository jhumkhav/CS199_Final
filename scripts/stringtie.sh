#!/bin/bash
#$ -N stringtie
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load stringtie

stringtie -p 8 -G genes.gtf -o C1_R1.gtf -l C1_R1 C1_R1.bam
stringtie -p 8 -G genes.gtf -o C1_R2.gtf -l C1_R2 C1_R2.bam
stringtie -p 8 -G genes.gtf -o C1_R3.gtf -l C1_R3 C1_R3.bam
stringtie -p 8 -G genes.gtf -o C2_R1.gtf -l C2_R1 C2_R1.bam
stringtie -p 8 -G genes.gtf -o C2_R2.gtf -l C2_R2 C2_R2.bam
stringtie -p 8 -G genes.gtf -o C2_R3.gtf -l C2_R3 C2_R3.bam
