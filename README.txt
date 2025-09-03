WES & RNA-seq Variant Analysis Workflow
=======================================

Introduction
------------

This workflow provides an end-to-end framework for analyzing high-throughput sequencing data using Snakemake, a scalable workflow management system that ensures reproducibility and efficient use of computational resources. The pipeline is designed primarily for whole-exome sequencing (WES) and RNA-seq datasets, enabling comprehensive exploration of genomic and transcriptomic variation. Key steps include raw data quality control, alignment, and variant calling from DNA and RNA reads, followed by integrated prioritization of variants across DNA and RNA to highlight potentially functional alterations. In addition, the pipeline supports allele-specific expression (ASE) analysis, which leverages RNA-seq to assess the transcriptional impact of genetic variants. By combining genomic and transcriptomic evidence, this workflow provides a systematic approach for identifying biologically relevant variants with potential regulatory or functional consequences.

* * *

Features
--------

*   Automated pipeline management using Snakemake
*   Support for WES and RNA-seq data
*   Quality control (FastQC, MultiQC)
*   Preprocessing (adapter trimming, filtering)
*   Alignment using BWA-MEM (DNA) and STAR (RNA)
*   Variant calling with GATK / bcftools
*   Variant prioritization across DNA and RNA
*   Allele-specific expression (ASE) analysis
*   Scalable execution on HPC or cloud environments

## Installation
### Requirements
- Snakemake (>=7.0)
- Conda (>=4.8)
- Python (>=3.8)

# Install Miniconda (recommended lightweight Conda)
#Download Miniconda installer for Linux (Python 3.8+)
# For other OS: see https://docs.conda.io/en/latest/miniconda.html
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts, then activate conda
source ~/.bashrc

# Update Conda to the latest version (>=4.8)
conda update -n base -c defaults conda

# Create a new environment with Python >=3.8
conda create -n snakemake python=3.9 -y
conda activate snakemake

# Install Snakemake (>=7.0)
conda install -c conda-forge -c bioconda snakemake -y

# Check versions
python --version
conda --version
snakemake --version


### Setup
Clone the repository:

```bash
git clone https://github.com/your-repo/ExomeRNA-VariantPipeline.git
cd ExomeRNA-VariantPipeline
```


