Whole Genome Sequencing (WGS) Variant Calling Pipeline

Overview

This repository contains a Bash-based pipeline for variant calling from Whole Genome Sequencing (WGS) FASTQ files. The workflow performs alignment, duplicate marking, base quality score recalibration (BQSR), and variant calling using widely used bioinformatics tools such as BWA, SAMtools, and GATK.

The pipeline aligns paired-end sequencing reads to the human reference genome (hg38) and generates a VCF file containing detected variants.

Pipeline Workflow

The pipeline follows the standard best practices for variant discovery.

1. Download Input Data

Patient paired-end FASTQ files

Reference genome and required annotation resources

Example input files: 
Sample_1.fq.gz/ Sample_1.fastq.gz
Sample_2.fq.gz/ Sample_2.fastq.gz

Reference Genome Preparation
Download Reference Genome

The human reference genome is downloaded from the European Bioinformatics Institute database.

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
Decompress and Rename
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa hg38.fa
