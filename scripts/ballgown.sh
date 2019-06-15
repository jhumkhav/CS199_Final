#!/bin/bash
#$ -N ballgown
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load stringtie

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R1/C1_R1.gtf C1_R1.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R2/C1_R2.gtf C1_R2.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R3/C1_R3.gtf C1_R3.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R1/C2_R1.gtf C2_R1.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R2/C2_R2.gtf C2_R2.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R3/C2_R3.gtf C2_R3.bam
