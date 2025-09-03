
configfile: "config/config.yaml"

import pandas as pd

# Load sample table and extract Run IDs
SAMPLES = pd.read_csv(config["sample_table"])
SRA_IDS = SAMPLES["Run"].tolist()
SAMPLE_SOURCE = dict(zip(SAMPLES["Run"], SAMPLES["LibrarySource"].str.lower()))
GENOMIC_IDS = SAMPLES.loc[SAMPLES["LibrarySource"].str.lower() == "genomic", "Run"].tolist()
TRANSCRIPTOMIC_IDS = SAMPLES.loc[SAMPLES["LibrarySource"].str.lower() == "transcriptomic", "Run"].tolist()


include: "rules/qc.smk"
include: "rules/get_data.smk"

rule all:
    input:
        # Downloaded data from SRA
        expand("data/raw_reads/{source}_reads/{id}_1.fastq.gz",
               id=SRA_IDS, source=[SAMPLE_SOURCE[i] for i in SRA_IDS]),
        expand("data/raw_reads/{source}_reads/{id}_2.fastq.gz",
               id=SRA_IDS, source=[SAMPLE_SOURCE[i] for i in SRA_IDS]),

        # FastQC on raw SRA reads (with source subfolders)
        expand("output/qc/fastqc/{source}/{id}_fastqc.html",
               id=SRA_IDS, source=[SAMPLE_SOURCE[i] for i in SRA_IDS]),
        expand("output/qc/fastqc/{source}/{id}_fastqc.zip",
               id=SRA_IDS, source=[SAMPLE_SOURCE[i] for i in SRA_IDS]),

        # MultiQC (separate for genomic + transcriptomic)
        "output/qc/multiqc/genomic/multiqc_report.html",
        "output/qc/multiqc/genomic/failed_samples.txt",
        "output/qc/multiqc/transcriptomic/multiqc_report.html",
        "output/qc/multiqc/transcriptomic/failed_samples.txt"

