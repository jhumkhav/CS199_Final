# Differential Gene and Transcript Expression Analysis

## Background
RNA-seq is used in order to measure and compare gene expression across different situations or conditions. This project aims to compare two different protocols that have been used to study differential gene expression. The first protocol uses TopHat and Cufflinks, while the second uses HISAT, StringTie, and Ballgown.

## Experimental Design
### Downloading Data
**TopHat and Cufflinks**
Download and unpack the Fruit Fly iGenome data.
'''
wget ftp://igenome:G3nom3s4u@ussdftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5.25/Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
tar -zxvf Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
'''
Create links to files that will be needed in the protocol - genes.gtf, genome.fa, and the Bowtie Index.
'''
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/WholeGenomeFasta/genome.fa

ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/WholeGenomeFasta/genome.fa
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.1.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.2.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.3.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.4.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.1.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.2.bt2
'''

Download and unpack the sequencing data, which is available through the Gene Expression Omnibus accession GSE32038.
'''
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038%5Fsimulated%5Ffastq%5Ffiles%2Etar%2Egz
tar -zxvf GSE32038_simulated_fastq_files.tar.gz
'''

**HISAT, StringTie, and Ballgown**
Unpack the data available through the following link:
ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
'''
tar xvzf chrX_data.tar.gz
'''

### Running the Analysis
**TopHat and Cufflinks**
Load the TopHat module and map the reads for each of the samples to the reference genome.
'''
module load tophat

tophat -p 8 -G genes.gtf -o C1_R1_thout genome GSM794483_C1_R1_1.fq GSM794483_C1_R1_2.fq
tophat -p 8 -G genes.gtf -o C1_R2_thout genome GSM794484_C1_R2_1.fq GSM794484_C1_R2_2.fq
tophat -p 8 -G genes.gtf -o C1_R3_thout genome GSM794485_C1_R3_1.fq GSM794485_C1_R3_2.fq
tophat -p 8 -G genes.gtf -o C2_R1_thout genome GSM794486_C2_R1_1.fq GSM794486_C2_R1_2.fq
tophat -p 8 -G genes.gtf -o C2_R2_thout genome GSM794487_C2_R2_1.fq GSM794487_C2_R2_2.fq
tophat -p 8 -G genes.gtf -o C2_R3_thout genome GSM794488_C2_R3_1.fq GSM794488_C2_R3_2.fq
'''

Load the Cufflinks module and assemble the transcripts for each sample.
'''
module load cufflinks

cufflinks -p 8 -o C1_R1_clout C1_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R2_clout C1_R2_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R3_clout C1_R3_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R1_clout C2_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R2_clout C2_R2_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R3_clout C2_R3_thout/accepted_hits.bam
'''

Create a file - assemblies.txt - that contains the assembly file for all of the samples, using vim. 
'''
vim assemblies
#Include the following text in the file:
#./C1_R1_clout/transcripts.gtf
#./C1_R2_clout/transcripts.gtf
#./C1_R3_clout/transcripts.gtf
#./C2_R1_clout/transcripts.gtf
#./C2_R2_clout/transcripts.gtf
#./C2_R3_clout/transcripts.gtf
'''

Load Cufflinks and create a single annotation that merges all of the assemblies.
'''
module load cufflinks

cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies.txt
'''

Use the merged assembly with the BAM files to identify differentially expressed genes and trancripts.
'''
module load cufflinks

cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794483_C1_R1_thout/accepted_hits.bam, GSM794486_C2_R1_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794484_C1_R2_thout/accepted_hits.bam, GSM794487_C2_R2_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794485_C1_R3_thout/accepted_hits.bam, GSM794488_C2_R3_thout/accepted_hits.bam
'''
