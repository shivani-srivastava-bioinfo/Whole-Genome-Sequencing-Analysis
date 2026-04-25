# Whole Genome Sequencing Analysis

A complete, reproducible germline variant calling pipeline following **GATK4 Best Practices** — from raw paired-end FASTQ files to variant-called VCF, with every step documented and tested.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Reference & Resource Files](#reference--resource-files)
- [Usage](#usage)
- [Pipeline Steps Detail](#pipeline-steps-detail)
  - [Step 1: Reference Preparation](#step-1-reference-preparation)
  - [Step 2: Read Alignment — BWA-MEM](#step-2-read-alignment--bwa-mem)
  - [Step 3: Mark Duplicates & Sort — GATK4](#step-3-mark-duplicates--sort--gatk4)
  - [Step 4: Base Quality Score Recalibration (BQSR)](#step-4-base-quality-score-recalibration-bqsr)
  - [Step 5: Variant Calling — HaplotypeCaller](#step-5-variant-calling--haplotypecaller)
- [Output Files](#output-files)
- [Notes & Recommendations](#notes--recommendations)
- [Repository Structure](#repository-structure)
- [Tool Versions](#tool-versions)
- [References](#references)

---

## Overview

This pipeline performs end-to-end germline short variant discovery (SNPs and INDELs) from Illumina paired-end whole-genome sequencing data. It covers:

- Reference genome indexing and preparation
- Paired-end read alignment to GRCh38/hg38 using **BWA-MEM**
- Duplicate marking and coordinate sorting using **GATK MarkDuplicatesSpark**
- Base quality recalibration (**BQSR**) using known variant sites from the GATK Resource Bundle
- Variant calling using **GATK HaplotypeCaller**

The pipeline is designed for single-sample germline analysis and produces a raw VCF ready for downstream filtering and annotation.

---

## Pipeline Architecture

```
Raw FASTQ (R1 + R2)
        │
        ▼
  BWA-MEM Alignment  ──────────────────────────► urgent.paired.sam
        │
        ▼
  MarkDuplicatesSpark  ────────────────────────► patient_sorted_dedup_reads.bam
  (sort + deduplicate)
        │
        ▼
  BaseRecalibrator  ───────────────────────────► recal_data.table
  (build BQSR model using dbSNP138)
        │
        ▼
  ApplyBQSR  ──────────────────────────────────► p_1_sorted_dedup_bqsr_reads.bam
        │
        ▼
  HaplotypeCaller  ────────────────────────────► raw_variants.vcf
```

---

## Requirements

Ensure the following tools are installed and available in your `$PATH` before running the pipeline:

| Tool | Version Recommended | Purpose |
|------|---------------------|---------|
| [BWA](https://github.com/lh3/bwa) | ≥ 0.7.17 | Paired-end read alignment |
| [GATK4](https://gatk.broadinstitute.org/) | ≥ 4.3.0 | Duplicate marking, BQSR, variant calling |
| [SAMtools](http://www.htslib.org/) | ≥ 1.15 | Reference FASTA indexing |
| [wget](https://www.gnu.org/software/wget/) | Any | Downloading resource bundle files |
| Java | 8 or 11 | Required by GATK4 |

> **System Requirements:** ≥ 16 GB RAM recommended. ≥ 200 GB free disk space for intermediate files (SAM, BAM, BWA index).

---

## Input Files

Place the following paired-end FASTQ files in the same directory as the script:

| File | Description |
|------|-------------|
| `Sample_1.fq.gz/ Sample_1.fastq.gz` | Read 1 — Forward reads |
| `Sample_2.fq.gz/ Sample_2.fastq.gz` | Read 2 — Reverse reads |

> The sample ID `Sample` is embedded into the BAM read group tag. Update the `-R` flag in the BWA-MEM command if using a different sample.

---

## Reference & Resource Files

### Reference Genome (hg38 / GRCh38)

The reference genome must be downloaded and decompressed manually before running the script:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa hg38.fa
```

The pipeline will automatically generate all required index and dictionary files:

| File | Tool | Purpose |
|------|------|---------|
| `hg38.fa.fai` | SAMtools faidx | FASTA index for random access |
| `hg38.dict` | GATK CreateSequenceDictionary | Sequence dictionary for GATK tools |
| `hg38.fa.{amb,ann,bwt,pac,sa}` | BWA index | Alignment index files |

### Known Sites — BQSR Resource Bundle

Downloaded automatically by the script from the GATK Resource Bundle (Google Cloud):

```
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

These files are used by `BaseRecalibrator` to distinguish true variants from sequencing errors during base quality score recalibration.

---

## Usage

Once all input and reference files are ready, run the script:

```bash
bash variant_calling_pipeline.sh
```

To capture all logs:

```bash
bash variant_calling_pipeline.sh 2>&1 | tee pipeline_run.log
```

---

## Pipeline Steps Detail

### Step 1: Reference Preparation

```bash
# Decompress and rename reference
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa hg38.fa

# Create .fai index (required by GATK)
samtools faidx hg38.fa

# Create .dict sequence dictionary (required by GATK)
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

# Download known variant sites for BQSR
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

**Why:** GATK tools require both a `.fai` index and `.dict` dictionary alongside the reference FASTA. Without these, HaplotypeCaller and other GATK tools will fail to start.

---

### Step 2: Read Alignment — BWA-MEM

```bash
# Index reference for BWA
bwa index hg38.fa

# Align paired-end reads
bwa mem -t 40 \
  -R "@RG\tID:Sample\tPL:ILLUMINA\tSM:Sample" \
  hg38.fa \
  Sample_1.fq.gz \
  Sample_2.fq.gz \
  > urgent.paired.sam
```

**Key flags:**

| Flag | Meaning |
|------|---------|
| `-t 40` | Use 40 CPU threads |
| `-R` | Read group tag — required by GATK downstream tools |

**Read Group fields:**

| Field | Value | Meaning |
|-------|-------|---------|
| `ID` | Sample | Unique run/lane identifier |
| `PL` | ILLUMINA | Sequencing platform |
| `SM` | Sample | Sample name |

> **Why Read Groups matter:** GATK HaplotypeCaller will throw an error if read groups are missing or incomplete. Always include `ID`, `PL`, and `SM` at minimum.

---

### Step 3: Mark Duplicates & Sort — GATK4

```bash
gatk MarkDuplicatesSpark \
  -I urgent.paired.sam \
  -O patient_sorted_dedup_reads.bam
```

**What it does:**
- Identifies and flags PCR duplicate reads (reads that originated from the same DNA template)
- Simultaneously sorts the output BAM by coordinate
- Uses Spark for parallelized execution

> **Why remove duplicates:** PCR duplicates inflate variant allele frequencies and can introduce false-positive variant calls. Marking them tells GATK to ignore duplicates during variant calling.

---

### Step 4: Base Quality Score Recalibration (BQSR)

BQSR is a two-pass process that corrects systematic errors in the base quality scores assigned by the sequencing instrument.

**Pass 1 — Build the recalibration model:**

```bash
gatk BaseRecalibrator \
  -I patient_sorted_dedup_reads.bam \
  -R hg38.fa \
  --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
  -O recal_data.table
```

**Pass 2 — Apply the model:**

```bash
gatk ApplyBQSR \
  -I patient_sorted_dedup_reads.bam \
  -R hg38.fa \
  --bqsr-recal-file recal_data.table \
  -O p_1_sorted_dedup_bqsr_reads.bam
```

> **Why BQSR:** Sequencing machines systematically over- or under-estimate base call quality depending on context (e.g., GC content, position in read, neighboring bases). BQSR detects and corrects these biases by comparing observed vs. expected error rates at known variant positions, resulting in more accurate variant calls downstream.

---

### Step 5: Variant Calling — HaplotypeCaller

```bash
gatk HaplotypeCaller \
  -R hg38.fa \
  -I p_1_sorted_dedup_bqsr_reads.bam \
  -O raw_variants.vcf
```

**What it does:**
- Calls SNPs and INDELs simultaneously using local de novo assembly
- Re-assembles active regions of the genome around candidate variant sites
- Applies a pair-HMM model to compute genotype likelihoods
- Outputs a raw VCF containing all called variant sites

> **For multi-sample cohort analysis:** Run HaplotypeCaller in GVCF mode (`-ERC GVCF`) per sample, then use `GenomicsDBImport` + `GenotypeGVCFs` for joint genotyping.

---

## Output Files

| File | Description |
|------|-------------|
| `urgent.paired.sam` | Raw BWA-MEM alignment (SAM format) |
| `patient_sorted_dedup_reads.bam` | Coordinate-sorted, duplicate-marked BAM |
| `recal_data.table` | BQSR recalibration model (text table) |
| `p_1_sorted_dedup_bqsr_reads.bam` | Final analysis-ready BAM (post-BQSR) |
| `raw_variants.vcf` | Raw germline variant calls — SNPs + INDELs |

---

## Notes & Recommendations

**Thread tuning:** Adjust `-t 40` in BWA-MEM to match available CPU cores. For `MarkDuplicatesSpark`, add `--spark-master local[4]` to parallelize across 40 cores.

**Disk space:** Ensure at least **200 GB** of free disk. The uncompressed SAM file alone can exceed 50–100 GB for WGS at 30× coverage.

---

## Repository Structure

```
variant-calling-pipeline/
├── variant_calling_pipeline.sh    ← Main pipeline script
├── README.md                      ← This file
└── logs/                          ← Log output directory (create before running)
```

---

## Tool Versions

| Tool | Version | Step |
|------|---------|------|
| BWA | ≥ 0.7.17 | Step 2 — Alignment |
| GATK4 | ≥ 4.3.0 | Steps 3, 4, 5 |
| SAMtools | ≥ 1.15 | Step 1 — Reference indexing |
| Java | 8 or 11 | GATK4 dependency |

---

## References

- **GATK Best Practices — Germline Short Variant Discovery:**
  https://gatk.broadinstitute.org/hc/en-us/articles/360035535932

- **BWA-MEM:**
  Li H. & Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Aligner. *Bioinformatics*, 25(14), 1754–1760.

- **GATK — HaplotypeCaller:**
  McKenna A. et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297–1303.

- **GATK Resource Bundle (hg38):**
  https://gatk.broadinstitute.org/hc/en-us/articles/360035890811

- **GENCODE Human Reference Genome Release 38:**
  https://www.gencodegenes.org/human/

- **SAMtools:**
  Danecek P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.

---

*Pipeline developed for single-sample germline WGS analysis following GATK4 Best Practices. All commands verified for GRCh38/hg38.*
