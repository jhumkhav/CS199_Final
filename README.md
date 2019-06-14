# Differential Gene and Transcript Expression Analysis

## Background
RNA-seq is used in order to measure and compare gene expression across different situations or conditions. This project aims to compare two different protocols that have been used to study differential gene expression. The first protocol uses TopHat and Cufflinks, while the second uses HISAT, StringTie, and Ballgown.

## Experimental Design
### TopHat and Cufflinks Protocol
**Download Data**

Download and unpack the Fruit Fly iGenome data.
```
wget ftp://igenome:G3nom3s4u@ussdftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5.25/Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
tar -zxvf Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
```

Create links to files that will be needed in the protocol - genes.gtf, genome.fa, and the Bowtie Index.
```
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/WholeGenomeFasta/genome.fa

ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/WholeGenomeFasta/genome.fa
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.1.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.2.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.3.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.4.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.1.bt2
ln -s ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.2.bt2
```

Download and unpack the sequencing data, which is available through the Gene Expression Omnibus accession GSE32038.
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038%5Fsimulated%5Ffastq%5Ffiles%2Etar%2Egz
tar -zxvf GSE32038_simulated_fastq_files.tar.gz
```

**Analysis**

Load the TopHat module and map the reads for each of the samples to the reference genome.
```
module load tophat

tophat -p 8 -G genes.gtf -o C1_R1_thout genome GSM794483_C1_R1_1.fq GSM794483_C1_R1_2.fq
tophat -p 8 -G genes.gtf -o C1_R2_thout genome GSM794484_C1_R2_1.fq GSM794484_C1_R2_2.fq
tophat -p 8 -G genes.gtf -o C1_R3_thout genome GSM794485_C1_R3_1.fq GSM794485_C1_R3_2.fq
tophat -p 8 -G genes.gtf -o C2_R1_thout genome GSM794486_C2_R1_1.fq GSM794486_C2_R1_2.fq
tophat -p 8 -G genes.gtf -o C2_R2_thout genome GSM794487_C2_R2_1.fq GSM794487_C2_R2_2.fq
tophat -p 8 -G genes.gtf -o C2_R3_thout genome GSM794488_C2_R3_1.fq GSM794488_C2_R3_2.fq
```

Load the Cufflinks module and assemble the transcripts for each sample.
```
module load cufflinks

cufflinks -p 8 -o C1_R1_clout C1_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R2_clout C1_R2_thout/accepted_hits.bam
cufflinks -p 8 -o C1_R3_clout C1_R3_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R1_clout C2_R1_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R2_clout C2_R2_thout/accepted_hits.bam
cufflinks -p 8 -o C2_R3_clout C2_R3_thout/accepted_hits.bam
```

Create a file - assemblies.txt - that contains the assembly file for all of the samples, using vim. 
```
vim assemblies
#Include the following text in the file:
#./C1_R1_clout/transcripts.gtf
#./C1_R2_clout/transcripts.gtf
#./C1_R3_clout/transcripts.gtf
#./C2_R1_clout/transcripts.gtf
#./C2_R2_clout/transcripts.gtf
#./C2_R3_clout/transcripts.gtf
```

Load Cufflinks and create a single annotation that merges all of the assemblies.
```
module load cufflinks

cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies.txt
```

Use the merged assembly with the BAM files to identify differentially expressed genes and trancripts.
```
module load cufflinks

cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794483_C1_R1_thout/accepted_hits.bam, GSM794486_C2_R1_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794484_C1_R2_thout/accepted_hits.bam, GSM794487_C2_R2_thout/accepted_hits.bam
cuffdiff -o diff_out -b genome.fa -L C1,C2 -u merged_asm/merged.gtf GSM794485_C1_R3_thout/accepted_hits.bam, GSM794488_C2_R3_thout/accepted_hits.bam
```

**Directory Tree**
Below is a tree of the directory for this protocol that contains all of the downloaded and created files.


### HISAT, StringTie, and Ballgown
**Downloading Data**

Unpack the data available through the following link:
ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
```
tar xvzf chrX_data.tar.gz
```
**Analysis**

Load the HISAT2 module and map the reads for each of the samples to the reference genome.
```
module load hisat2

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188104_chrX_1.fastq.gz -2 chrX_data/samples/ERR188104_chrX_2.fastq.gz -S ERR188104_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188234_chrX_1.fastq.gz -2 chrX_data/samples/ERR188234_chrX_2.fastq.gz -S ERR188234_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188245_chrX_1.fastq.gz -2 chrX_data/samples/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188257_chrX_1.fastq.gz -2 chrX_data/samples/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188273_chrX_1.fastq.gz -2 chrX_data/samples/ERR188273_chrX_2.fastq.gz -S ERR188273_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188337_chrX_1.fastq.gz -2 chrX_data/samples/ERR188337_chrX_2.fastq.gz -S ERR188337_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188383_chrX_1.fastq.gz -2 chrX_data/samples/ERR188383_chrX_2.fastq.gz -S ERR188383_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188401_chrX_1.fastq.gz -2 chrX_data/samples/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188454_chrX_1.fastq.gz -2 chrX_data/samples/ERR188454_chrX_2.fastq.gz -S ERR188454_chrX.sam
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR204916_chrX_1.fastq.gz -2 chrX_data/samples/ERR204916_chrX_2.fastq.gz -S ERR204916_chrX.sam
```

Load samtools module and convert SAM files to BAM files
```
module load samtools

samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam
samtools sort -@ 8 -o ERR188104_chrX.bam ERR188104_chrX.sam
samtools sort -@ 8 -o ERR188234_chrX.bam ERR188234_chrX.sam
samtools sort -@ 8 -o ERR188245_chrX.bam ERR188245_chrX.sam
samtools sort -@ 8 -o ERR188257_chrX.bam ERR188257_chrX.sam
samtools sort -@ 8 -o ERR188273_chrX.bam ERR188273_chrX.sam
samtools sort -@ 8 -o ERR188337_chrX.bam ERR188337_chrX.sam
samtools sort -@ 8 -o ERR188383_chrX.bam ERR188383_chrX.sam
samtools sort -@ 8 -o ERR188401_chrX.bam ERR188401_chrX.sam
samtools sort -@ 8 -o ERR188428_chrX.bam ERR188428_chrX.sam
samtools sort -@ 8 -o ERR188454_chrX.bam ERR188454_chrX.sam
samtools sort -@ 8 -o ERR204916_chrX.bam ERR204916_chrX.sam
```

Load StringTie and assemble the transcripts for each file.
```
module load stringtie

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 ERR188044_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188104_chrX.gtf -l ERR188104 ERR188104_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188234_chrX.gtf -l ERR188234 ERR188234_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188245_chrX.gtf -l ERR188245 ERR188245_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188257_chrX.gtf -l ERR188257 ERR188257_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188273_chrX.gtf -l ERR188273 ERR188273_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188337_chrX.gtf -l ERR188337 ERR188337_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188383_chrX.gtf -l ERR188383 ERR188383_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188401_chrX.gtf -l ERR188401 ERR188401_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188428_chrX.gtf -l ERR188428 ERR188428_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188454_chrX.gtf -l ERR188454 ERR188454_chrX.bam
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR204916_chrX.gtf -l ERR204916 ERR204916_chrX.bam
```

Merge all of the trancripts using StringTie module.
```
module load stringtie

stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
```

Load StringeTie module and determine abundances of transcripts and prepare data for Ballgown.
```
module load stringtie

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188104/ERR188104_chrX.gtf ERR188104_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188234/ERR188234_chrX.gtf ERR188234_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188245/ERR188245_chrX.gtf ERR188245_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188257/ERR188257_chrX.gtf ERR188257_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188273/ERR188273_chrX.gtf ERR188273_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188337/ERR188337_chrX.gtf ERR188337_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188383/ERR188383_chrX.gtf ERR188383_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188401/ERR188401_chrX.gtf ERR188401_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188428/ERR188428_chrX.gtf ERR188428_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188454/ERR188454_chrX.gtf ERR188454_chrX.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR204916/ERR204916_chrX.gtf ERR204916_chrX.bam
```

**Directory Tree**
Below is a tree of the directory for this protocol that contains all of the downloaded and created files.
```
|-- ballgown
|   |-- ERR188044
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188044_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188104
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188104_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188234
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188234_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188245
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188245_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188257
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188257_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188273
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188273_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188337
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188337_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188383
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188383_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188401
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188401_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188428
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188428_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- ERR188454
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- ERR188454_chrX.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   `-- ERR204916
|       |-- e2t.ctab
|       |-- e_data.ctab
|       |-- ERR204916_chrX.gtf
|       |-- i2t.ctab
|       |-- i_data.ctab
|       `-- t_data.ctab
|-- ballgown.sh
|-- bdgp6
|   |-- genome.1.ht2
|   |-- genome.2.ht2
|   |-- genome.3.ht2
|   |-- genome.4.ht2
|   |-- genome.5.ht2
|   |-- genome.6.ht2
|   |-- genome.7.ht2
|   |-- genome.8.ht2
|   `-- make_bdgp6.sh
|-- bdgp6.tar.gz
|-- chrX_data
|   |-- genes
|   |   `-- chrX.gtf
|   |-- genome
|   |   `-- chrX.fa
|   |-- geuvadis_phenodata.csv
|   |-- indexes
|   |   |-- chrX_tran.1.ht2
|   |   |-- chrX_tran.2.ht2
|   |   |-- chrX_tran.3.ht2
|   |   |-- chrX_tran.4.ht2
|   |   |-- chrX_tran.5.ht2
|   |   |-- chrX_tran.6.ht2
|   |   |-- chrX_tran.7.ht2
|   |   `-- chrX_tran.8.ht2
|   |-- mergelist.txt
|   `-- samples
|       |-- ERR188044_chrX_1.fastq.gz
|       |-- ERR188044_chrX_2.fastq.gz
|       |-- ERR188104_chrX_1.fastq.gz
|       |-- ERR188104_chrX_2.fastq.gz
|       |-- ERR188234_chrX_1.fastq.gz
|       |-- ERR188234_chrX_2.fastq.gz
|       |-- ERR188245_chrX_1.fastq.gz
|       |-- ERR188245_chrX_2.fastq.gz
|       |-- ERR188257_chrX_1.fastq.gz
|       |-- ERR188257_chrX_2.fastq.gz
|       |-- ERR188273_chrX_1.fastq.gz
|       |-- ERR188273_chrX_2.fastq.gz
|       |-- ERR188337_chrX_1.fastq.gz
|       |-- ERR188337_chrX_2.fastq.gz
|       |-- ERR188383_chrX_1.fastq.gz
|       |-- ERR188383_chrX_2.fastq.gz
|       |-- ERR188401_chrX_1.fastq.gz
|       |-- ERR188401_chrX_2.fastq.gz
|       |-- ERR188428_chrX_1.fastq.gz
|       |-- ERR188428_chrX_2.fastq.gz
|       |-- ERR188454_chrX_1.fastq.gz
|       |-- ERR188454_chrX_2.fastq.gz
|       |-- ERR204916_chrX_1.fastq.gz
|       `-- ERR204916_chrX_2.fastq.gz
|-- ERR188044_chrX.bam
|-- ERR188044_chrX.gtf
|-- ERR188044_chrX.sam
|-- ERR188104_chrX.bam
|-- ERR188104_chrX.gtf
|-- ERR188104_chrX.sam
|-- ERR188234_chrX.bam
|-- ERR188234_chrX.gtf
|-- ERR188234_chrX.sam
|-- ERR188245_chrX.bam
|-- ERR188245_chrX.gtf
|-- ERR188245_chrX.sam
|-- ERR188257_chrX.bam
|-- ERR188257_chrX.gtf
|-- ERR188257_chrX.sam
|-- ERR188273_chrX.bam
|-- ERR188273_chrX.gtf
|-- ERR188273_chrX.sam
|-- ERR188337_chrX.bam
|-- ERR188337_chrX.gtf
|-- ERR188337_chrX.sam
|-- ERR188383_chrX.bam
|-- ERR188383_chrX.gtf
|-- ERR188383_chrX.sam
|-- ERR188401_chrX.bam
|-- ERR188401_chrX.gtf
|-- ERR188401_chrX.sam
|-- ERR188428_chrX.bam
|-- ERR188428_chrX.gtf
|-- ERR188428_chrX.sam
|-- ERR188454_chrX.bam
|-- ERR188454_chrX.gtf
|-- ERR188454_chrX.sam
|-- ERR204916_chrX.bam
|-- ERR204916_chrX.gtf
|-- ERR204916_chrX.sam
|-- hisat.sh
|-- merge.sh
|-- my_hisat.sh
|-- samtools.sh
|-- stringtie_merged.gtf
`-- stringtie.sh
