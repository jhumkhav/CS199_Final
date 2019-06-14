#!/bin/bash
#$ -N hisat
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load hisat2

hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794483_C1_R1_1.fq.gz -2 GSM794483_C1_R1_2.fq.gz -S C1_R1.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794484_C1_R2_1.fq.gz -2 GSM794484_C1_R2_2.fq.gz -S C1_R2.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794485_C1_R3_1.fq.gz -2 GSM794485_C1_R3_2.fq.gz -S C1_R3.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794486_C2_R1_1.fq.gz -2 GSM794486_C2_R1_2.fq.gz -S C2_R1.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794487_C2_R2_1.fq.gz -2 GSM794487_C2_R2_2.fq.gz -S C2_R2.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794488_C2_R3_1.fq.gz -2 GSM794488_C2_R3_2.fq.gz -S C2_R3.sam
