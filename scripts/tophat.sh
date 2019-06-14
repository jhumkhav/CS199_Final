
#!/bin/bash
#$ -N tophat
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 6

module load cufflinks
module load samtools
module load tophat
module load bowtie2

tophat -p 8 -G genes.gtf -o C1_R1_thout genome GSM794483_C1_R1_1.fq GSM794483_C1_R1_2.fq
tophat -p 8 -G genes.gtf -o C1_R2_thout genome GSM794484_C1_R2_1.fq GSM794484_C1_R2_2.fq
tophat -p 8 -G genes.gtf -o C1_R3_thout genome GSM794485_C1_R3_1.fq GSM794485_C1_R3_2.fq
tophat -p 8 -G genes.gtf -o C2_R1_thout genome GSM794486_C2_R1_1.fq GSM794486_C2_R1_2.fq
tophat -p 8 -G genes.gtf -o C2_R2_thout genome GSM794487_C2_R2_1.fq GSM794487_C2_R2_2.fq
tophat -p 8 -G genes.gtf -o C2_R3_thout genome GSM794488_C2_R3_1.fq GSM794488_C2_R3_2.fq

