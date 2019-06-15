#!/bin/bash
#$ -N merge
#$ -ckpt restart
#$ -q pub8i
#$ -pe openmp 1

module load stringtie

stringtie --merge -p 8 -G genes.gtf -o stringtie_merged.gtf mergelist.txt
