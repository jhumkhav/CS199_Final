# Differential Gene and Transcript Expression Analysis

## Background
RNA-seq is used in order to measure and compare gene expression across different situations or conditions. This project aims to compare two different protocols that have been used to study differential gene expression. The first protocol uses TopHat and Cufflinks, while the second uses HISAT, StringTie, and Ballgown.

## Experimental Design
**Create the Environments**

The necessary environments for this project are [TopHat](envs/tophat.yml), [Bowtie](envs/bowtie.yml), [Cufflinks](envs/cufflinks.yml), [HISAT](envs/hisat.yml), [StringTie](envs/stringtie.yml), and [Samtools](envs/samtools.yml).

Once the Miniconda3 was [installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), the environments listed above were created following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands).

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

**Data Visualization**

Loading R and the CummerBund package, create a database from the data in the Cuffdiff output.
```
> library(cummeRbund)
> cuff_data <- readCufflinks('diff_out')
```

Create a plot of the distribution of the expression levels for the samples.
```
> csDensity(genes(cuff_data))
```

Create a scatter plot that compares the expression of the genes in the different conditions.
```
> csScatter(genes(cuff_data), 'C1', 'C2')
```

**Directory Tree**

Below is a tree of the directory for this protocol that contains all of the downloaded and created files.
```
|-- assemblies.txt
|-- C1_R1_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C1_R1_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.merge_bam.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- reports.samtools_sort.log1
|   |   |-- reports.samtools_sort.log2
|   |   |-- reports.samtools_sort.log3
|   |   |-- reports.samtools_sort.log4
|   |   |-- reports.samtools_sort.log5
|   |   |-- reports.samtools_sort.log6
|   |   |-- reports.samtools_sort.log7
|   |   |-- run.log
|   |   `-- tophat.log
|   |-- prep_reads.info
|   `-- unmapped.bam
|-- C1_R2_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C1_R2_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.merge_bam.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- reports.samtools_sort.log1
|   |   |-- reports.samtools_sort.log2
|   |   |-- reports.samtools_sort.log3
|   |   |-- reports.samtools_sort.log4
|   |   |-- reports.samtools_sort.log5
|   |   |-- reports.samtools_sort.log6
|   |   |-- reports.samtools_sort.log7
|   |   |-- run.log
|   |   `-- tophat.log
|   |-- prep_reads.info
|   `-- unmapped.bam
|-- C1_R3_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C1_R3_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bowtie_build.log
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- bowtie.right_kept_reads.m2g_um.log
|   |   |-- bowtie.right_kept_reads.m2g_um_unmapped.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- juncs_db.log
|   |   |-- long_spanning_reads.segs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- run.log
|   |   `-- tophat.log
|   `-- prep_reads.info
|-- C2_R1_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C2_R1_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bam_merge_um.log
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.merge_bam.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- reports.samtools_sort.log1
|   |   |-- reports.samtools_sort.log2
|   |   |-- reports.samtools_sort.log3
|   |   |-- reports.samtools_sort.log4
|   |   |-- reports.samtools_sort.log5
|   |   |-- reports.samtools_sort.log6
|   |   |-- reports.samtools_sort.log7
|   |   |-- run.log
|   |   `-- tophat.log
|   |-- prep_reads.info
|   `-- unmapped.bam
|-- C2_R2_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C2_R2_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.merge_bam.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- reports.samtools_sort.log1
|   |   |-- reports.samtools_sort.log2
|   |   |-- reports.samtools_sort.log3
|   |   |-- reports.samtools_sort.log4
|   |   |-- reports.samtools_sort.log5
|   |   |-- reports.samtools_sort.log6
|   |   |-- reports.samtools_sort.log7
|   |   |-- run.log
|   |   `-- tophat.log
|   |-- prep_reads.info
|   `-- unmapped.bam
|-- C2_R3_clout
|   |-- genes.fpkm_tracking
|   |-- isoforms.fpkm_tracking
|   |-- skipped.gtf
|   `-- transcripts.gtf
|-- C2_R3_thout
|   |-- accepted_hits.bam
|   |-- align_summary.txt
|   |-- deletions.bed
|   |-- insertions.bed
|   |-- junctions.bed
|   |-- logs
|   |   |-- bowtie.left_kept_reads.log
|   |   |-- bowtie.right_kept_reads.log
|   |   |-- g2f.err
|   |   |-- g2f.out
|   |   |-- gtf_juncs.log
|   |   |-- m2g_left_kept_reads.err
|   |   |-- m2g_left_kept_reads.out
|   |   |-- m2g_right_kept_reads.err
|   |   |-- m2g_right_kept_reads.out
|   |   |-- prep_reads.log
|   |   |-- reports.log
|   |   |-- reports.merge_bam.log
|   |   |-- reports.samtools_sort.log0
|   |   |-- reports.samtools_sort.log1
|   |   |-- reports.samtools_sort.log2
|   |   |-- reports.samtools_sort.log3
|   |   |-- reports.samtools_sort.log4
|   |   |-- reports.samtools_sort.log5
|   |   |-- reports.samtools_sort.log6
|   |   |-- reports.samtools_sort.log7
|   |   |-- run.log
|   |   `-- tophat.log
|   `-- prep_reads.info
|-- cuffdiff.e7128987
|-- cuffdiff.o7128987
|-- cuffdiff.sh
|-- cufflinks2.sh
|-- cufflinks.sh
|-- cuffmerge.sh
|-- diff_out
|   |-- bias_params.info
|   |-- cds.count_tracking
|   |-- cds.diff
|   |-- cds_exp.diff
|   |-- cds.fpkm_tracking
|   |-- cds.read_group_tracking
|   |-- gene_exp.diff
|   |-- genes.count_tracking
|   |-- genes.fpkm_tracking
|   |-- genes.read_group_tracking
|   |-- isoform_exp.diff
|   |-- isoforms.count_tracking
|   |-- isoforms.fpkm_tracking
|   |-- isoforms.read_group_tracking
|   |-- promoters.diff
|   |-- read_groups.info
|   |-- run.info
|   |-- splicing.diff
|   |-- tss_group_exp.diff
|   |-- tss_groups.count_tracking
|   |-- tss_groups.fpkm_tracking
|   |-- tss_groups.read_group_tracking
|   `-- var_model.info
|-- Drosophila_melanogaster
|   `-- Ensembl
|       `-- BDGP5.25
|           |-- Annotation
|           |   |-- Archives
|           |   |   |-- archive-2010-09-27-21-29-49
|           |   |   |   |-- ChromInfo.txt
|           |   |   |   |-- DATE.txt
|           |   |   |   |-- Drosophila_melanogaster.BDGP5.25.59.gtf
|           |   |   |   |-- README.txt
|           |   |   |   |-- refFlat.txt
|           |   |   |   |-- refFlat.txt.gz
|           |   |   |   |-- refGene.txt
|           |   |   |   |-- splice_sites_34
|           |   |   |   |   |-- exon_coords.txt
|           |   |   |   |   |-- splice_sites-34.fa
|           |   |   |   |   |-- splice_sites-34.fa.2bpb
|           |   |   |   |   |-- splice_sites-34.fa.idx
|           |   |   |   |   `-- splice_sites-34.fa.vld
|           |   |   |   `-- splice_sites_49
|           |   |   |       |-- exon_coords.txt
|           |   |   |       |-- splice_sites-49.fa
|           |   |   |       |-- splice_sites-49.fa.2bpb
|           |   |   |       |-- splice_sites-49.fa.idx
|           |   |   |       `-- splice_sites-49.fa.vld
|           |   |   |-- archive-2011-01-27-18-19-23
|           |   |   |   |-- Genes
|           |   |   |   |   |-- ChromInfo.txt
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   |-- refFlat.txt.gz.bak
|           |   |   |   |   `-- refGene.txt
|           |   |   |   |-- README.txt
|           |   |   |   |-- SmallRNA
|           |   |   |   |   |-- hairpin.fa
|           |   |   |   |   `-- mature.fa
|           |   |   |   `-- Variation
|           |   |   |-- archive-2011-08-30-21-38-10
|           |   |   |   |-- Genes
|           |   |   |   |   |-- ChromInfo.txt
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   |-- refFlat.txt.gz.bak
|           |   |   |   |   `-- refGene.txt
|           |   |   |   |-- README.txt
|           |   |   |   |-- SmallRNA
|           |   |   |   |   |-- mature.fa
|           |   |   |   |   `-- precursor.fa
|           |   |   |   `-- Variation
|           |   |   |-- archive-2012-03-09-03-11-08
|           |   |   |   |-- Genes
|           |   |   |   |   |-- ChromInfo.txt
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   |-- refFlat.txt.gz.bak
|           |   |   |   |   `-- refGene.txt
|           |   |   |   |-- README.txt
|           |   |   |   |-- SmallRNA
|           |   |   |   |   |-- mature.fa
|           |   |   |   |   `-- precursor.fa
|           |   |   |   `-- Variation
|           |   |   |-- archive-2013-03-06-11-13-59
|           |   |   |   |-- Genes
|           |   |   |   |   |-- ChromInfo.txt
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   |-- refFlat.txt.gz.bak
|           |   |   |   |   `-- refGene.txt
|           |   |   |   |-- README.txt
|           |   |   |   |-- SmallRNA
|           |   |   |   |   |-- mature.fa
|           |   |   |   |   `-- precursor.fa
|           |   |   |   `-- Variation
|           |   |   |-- archive-2014-05-23-16-02-55
|           |   |   |   |-- Genes
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   `-- refFlat.txt.gz.bak
|           |   |   |   |-- README.txt
|           |   |   |   `-- SmallRNA
|           |   |   |       |-- hairpin.fa
|           |   |   |       `-- mature.fa
|           |   |   |-- archive-2015-07-17-14-30-26
|           |   |   |   |-- Genes
|           |   |   |   |   |-- genes.gtf
|           |   |   |   |   |-- genes.gtf.bak
|           |   |   |   |   |-- refFlat.txt.gz
|           |   |   |   |   `-- refFlat.txt.gz.bak
|           |   |   |   |-- README.txt
|           |   |   |   `-- SmallRNA
|           |   |   |       |-- hairpin.fa
|           |   |   |       `-- mature.fa
|           |   |   `-- archive-current -> archive-2015-07-17-14-30-26
|           |   |-- Genes -> Archives/archive-current/Genes
|           |   |-- README.txt -> Archives/archive-current/README.txt
|           |   `-- SmallRNA -> Archives/archive-current/SmallRNA
|           `-- Sequence
|               |-- AbundantSequences
|               |   |-- adapter_contam1.fa
|               |   |-- M.fa -> ../Chromosomes/M.fa
|               |   |-- phix.fa
|               |   |-- polyA.fa
|               |   |-- polyC.fa
|               |   `-- ribosomal.fa
|               |-- Bowtie2Index
|               |   |-- genome.1.bt2
|               |   |-- genome.2.bt2
|               |   |-- genome.3.bt2
|               |   |-- genome.4.bt2
|               |   |-- genome.fa -> ../WholeGenomeFasta/genome.fa
|               |   |-- genome.rev.1.bt2
|               |   `-- genome.rev.2.bt2
|               |-- BowtieIndex
|               |   |-- genome.1.ebwt
|               |   |-- genome.2.ebwt
|               |   |-- genome.3.ebwt
|               |   |-- genome.4.ebwt
|               |   |-- genome.fa -> ../WholeGenomeFasta/genome.fa
|               |   |-- genome.rev.1.ebwt
|               |   `-- genome.rev.2.ebwt
|               |-- BWAIndex
|               |   |-- genome.fa -> version0.6.0/genome.fa
|               |   |-- genome.fa.amb -> version0.6.0/genome.fa.amb
|               |   |-- genome.fa.ann -> version0.6.0/genome.fa.ann
|               |   |-- genome.fa.bwt -> version0.6.0/genome.fa.bwt
|               |   |-- genome.fa.pac -> version0.6.0/genome.fa.pac
|               |   |-- genome.fa.sa -> version0.6.0/genome.fa.sa
|               |   |-- version0.5.x
|               |   |   |-- genome.fa -> ../../WholeGenomeFasta/genome.fa
|               |   |   |-- genome.fa.amb
|               |   |   |-- genome.fa.ann
|               |   |   |-- genome.fa.bwt
|               |   |   |-- genome.fa.pac
|               |   |   |-- genome.fa.rbwt
|               |   |   |-- genome.fa.rpac
|               |   |   |-- genome.fa.rsa
|               |   |   `-- genome.fa.sa
|               |   `-- version0.6.0
|               |       |-- genome.fa -> ../../WholeGenomeFasta/genome.fa
|               |       |-- genome.fa.amb
|               |       |-- genome.fa.ann
|               |       |-- genome.fa.bwt
|               |       |-- genome.fa.pac
|               |       `-- genome.fa.sa
|               |-- Chromosomes
|               |   |-- 2L.fa
|               |   |-- 2R.fa
|               |   |-- 3L.fa
|               |   |-- 3R.fa
|               |   |-- 4.fa
|               |   |-- M.fa
|               |   `-- X.fa
|               `-- WholeGenomeFasta
|                   |-- genome.dict
|                   |-- genome.fa
|                   |-- genome.fa.fai
|                   `-- GenomeSize.xml
|-- genes.gtf -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf
|-- genome.1.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.1.bt2
|-- genome.2.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.2.bt2
|-- genome.3.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.3.bt2
|-- genome.4.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.4.bt2
|-- genome.fa -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/WholeGenomeFasta/genome.fa
|-- genome.rev.1.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.1.bt2
|-- genome.rev.2.bt2 -> ./Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome.rev.2.bt2
|-- GSM794483_C1_R1_1.fq
|-- GSM794483_C1_R1_2.fq
|-- GSM794484_C1_R2_1.fq
|-- GSM794484_C1_R2_2.fq
|-- GSM794485_C1_R3_1.fq
|-- GSM794485_C1_R3_2.fq
|-- GSM794486_C2_R1_1.fq
|-- GSM794486_C2_R1_2.fq
|-- GSM794487_C2_R2_1.fq
|-- GSM794487_C2_R2_2.fq
|-- GSM794488_C2_R3_1.fq
|-- GSM794488_C2_R3_2.fq
|-- merged_asm
|   |-- logs
|   |   `-- run.log
|   `-- merged.gtf
`-- tophat.sh
```

### HISAT, StringTie, and Ballgown
**Downloading Data**

Unpack the same data that was obtained in the previous protocol.
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038%5Fsimulated%5Ffastq%5Ffiles%2Etar%2Egz
tar -zxvf GSE32038_simulated_fastq_files.tar.gz
```
**Analysis**

Load the HISAT2 module and map the reads for each of the samples to the reference genome.
```
module load hisat2

hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794483_C1_R1_1.fq.gz -2 GSM794483_C1_R1_2.fq.gz -S C1_R1.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794484_C1_R2_1.fq.gz -2 GSM794484_C1_R2_2.fq.gz -S C1_R2.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794485_C1_R3_1.fq.gz -2 GSM794485_C1_R3_2.fq.gz -S C1_R3.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794486_C2_R1_1.fq.gz -2 GSM794486_C2_R1_2.fq.gz -S C2_R1.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794487_C2_R2_1.fq.gz -2 GSM794487_C2_R2_2.fq.gz -S C2_R2.sam
hisat2 -p 8 --dta -x bdgp6/genome -1 GSM794488_C2_R3_1.fq.gz -2 GSM794488_C2_R3_2.fq.gz -S C2_R3.sam
```

Load samtools module and convert SAM files to BAM files
```
module load samtools

samtools sort -@ 8 -o C1_R1.bam C1_R1.sam
samtools sort -@ 8 -o C1_R2.bam C1_R2.sam
samtools sort -@ 8 -o C1_R3.bam C1_R3.sam
samtools sort -@ 8 -o C2_R1.bam C2_R1.sam
samtools sort -@ 8 -o C2_R2.bam C2_R2.sam
samtools sort -@ 8 -o C2_R3.bam C2_R3.sam
```

Load StringTie and assemble the transcripts for each file.
```
module load stringtie

stringtie -p 8 -G genes.gtf -o C1_R1.gtf -l C1_R1 C1_R1.bam
stringtie -p 8 -G genes.gtf -o C1_R2.gtf -l C1_R2 C1_R2.bam
stringtie -p 8 -G genes.gtf -o C1_R3.gtf -l C1_R3 C1_R3.bam
stringtie -p 8 -G genes.gtf -o C2_R1.gtf -l C2_R1 C2_R1.bam
stringtie -p 8 -G genes.gtf -o C2_R2.gtf -l C2_R2 C2_R2.bam
stringtie -p 8 -G genes.gtf -o C2_R3.gtf -l C2_R3 C2_R3.bam
```

Merge all of the trancripts using StringTie module.
```
module load stringtie

stringtie --merge -p 8 -G genes.gtf -o stringtie_merged.gtf mergelist.txt
```

Load StringeTie module and determine abundances of transcripts and prepare data for Ballgown.
```
module load stringtie

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R1/C1_R1.gtf C1_R1.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R2/C1_R2.gtf C1_R2.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C1_R3/C1_R3.gtf C1_R3.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R1/C2_R1.gtf C2_R1.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R2/C2_R2.gtf C2_R2.bam
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/C2_R3/C2_R3.gtf C2_R3.bam
```

**Directory Tree**

Below is a tree of the directory for this protocol that contains all of the downloaded and created files.

```
|-- ballgown
|   |-- C1_R1
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C1_R1.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- C1_R2
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C1_R2.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- C1_R3
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C1_R3.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- C2_R1
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C2_R1.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- C2_R2
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C2_R2.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
|   |-- C2_R3
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- C2_R3.gtf
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   `-- t_data.ctab
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
|-- C1_R1.bam
|-- C1_R1.gtf
|-- C1_R1.sam
|-- C1_R2.bam
|-- C1_R2.gtf
|-- C1_R2.sam
|-- C1_R3.bam
|-- C1_R3.gtf
|-- C1_R3.sam
|-- C2_R1.bam
|-- C2_R1.gtf
|-- C2_R1.sam
|-- C2_R2.bam
|-- C2_R2.gtf
|-- C2_R2.sam
|-- C2_R3.bam
|-- C2_R3.gtf
|-- C2_R3.sam
|-- genes.gtf -> /data/users/jhumkhav/rnaseq/Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf
|-- GSE32038_simulated_fastq_files.tar.gz
|-- GSM794483_C1_R1_1.fq.gz
|-- GSM794483_C1_R1_2.fq.gz
|-- GSM794484_C1_R2_1.fq.gz
|-- GSM794484_C1_R2_2.fq.gz
|-- GSM794485_C1_R3_1.fq.gz
|-- GSM794485_C1_R3_2.fq.gz
|-- GSM794486_C2_R1_1.fq.gz
|-- GSM794486_C2_R1_2.fq.gz
|-- GSM794487_C2_R2_1.fq.gz
|-- GSM794487_C2_R2_2.fq.gz
|-- GSM794488_C2_R3_1.fq.gz
|-- GSM794488_C2_R3_2.fq.gz
|-- merge.sh
|-- mergelist.txt
|-- myhisat.sh
|-- samtools.sh
|-- stringtie_merged.gtf
`-- stringtie.sh


