# agu79: Zebrafish Genome Analysis

## Introduction

This report describes the genome sequence results for the "agu79". The following pipeline includes all steps from FASTQ download to SNP annotation.

## Dependencies

### Download DNA sequence & GTF

- `seqkit` (2.1.0)

### SNPs calling

- `fastq-dump` (2.11.3)
- `fastp` (0.20.1)
- `bwa-mem` (0.7.17)
- `gatk` (4.4.0.0)
- `vcftools` (0.1.16)
- `samtools` (1.13)
- `snpEff` (4.3)

### Genomic characteristics analysis

- `R` (4.2.3)
  - `ggplot2` (3.4.2)
  - `openxlsx` (4.2.5.2)
  - `patchwork` (1.1.2)
  - `ggsignif` (0.6.4)
  - `dplyr` (1.1.2)
  - `ggbreak` (0.1.1)
  - `stringr` (1.5.0)
  - `UpSetR` (1.4.0)
  - `reshape2` (1.4.4)
  - `sets` (1.0-24)

## Download DNA sequence & GTF

Reference data is downloaded by the following command:

```bash
wget https://ftp.ensembl.org/pub/release-101/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
```

For fa file as reference, only Chr 1\~25 and MT are used:

```bash
for i in `cat chr_mt.txt`
do
  seqkit grep -p ${i} Danio_rerio.GRCz11.dna.primary_assembly.fa >> Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa
done
```

To download GTF:

```bash
wget https://ftp.ensembl.org/pub/release-109/gff3/danio_rerio/Danio_rerio.GRCz11.109.gff3.gz
```

## Directory Structure

```
agu79/
├── bam/                    # Processed BAM files
├── chr.txt                 # Chromosomes (1–25)
├── chr_mt.txt              # Chromosomes (1–25, MT)
├── direc2.txt              # List of required directories
├── SRA.txt                 # SRA accession IDs
├── rawdata/
│   └── test.txt            # List of sample IDs for trimming
├── bam/                    # Processed BAM files
├── depth/                  # Coverage depth outputs
├── rawdata/                # FASTQ files
├── sam/                    # SAM alignment outputs
├── trimed/                 # Trimmed FASTQ files
├── vcf/
│   ├── 01/ ~ 06/           # Per-sample variant files
│   ├── HIGH/               # SnpEff annotated HIGH-impact variants
│   ├── MOD/                # SnpEff annotated MODERATE-impact variants
│   └── fin_vcf/            # Merged final VCFs
├── snpscall.sh             # Main pipeline script
├── vcf_analysis.R          # R script for downstream analysis
└── Danio_rerio.*           # Reference genome & GTF files
```

## BioProject

- **Project ID**: PRJNA1180756 (agu79)

### SRA Samples

| Sample name            | Accession   | Status   |
| ---------------------- | ----------- | -------- |
| agu79\_heterozygous\_1 | SRR31189664 | released |
| agu79\_homozygous\_1   | SRR31189661 | released |
| agu79\_homozygous\_2   | SRR31189662 | released |
| agu79\_homozygous\_3   | SRR31189663 | released |
| agu79\_homozygous\_4   | SRR31189660 | released |
| agu79\_homozygous\_5   | SRR31189659 | released |

## How to Run

Run the entire SNP pipeline with:

```bash
bash snpscall.sh
```

Run the downstream R analysis with:

```bash
Rscript vcf_analysis.R
```

---

For any questions, please contact the maintainer.

