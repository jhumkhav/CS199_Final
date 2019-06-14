#!/bin/bash
#$ -N cuffmerge
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load cufflinks

cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies.txt
