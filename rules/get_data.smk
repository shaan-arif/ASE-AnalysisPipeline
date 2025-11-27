# ==========================
# Rule: Download SRA Runs from BioProject PRJNA675090
# Purpose: Download raw SRA sequencing reads (FASTQ) for each sample
# from a given BioProject (e.g., PRJNA675090) using the SRA Run IDs.
# The reads are split into paired-end files and compressed with gzip.
# ==========================

rule get_data:
    output:
        "data/raw_reads/{source}_reads/{id}_1.fastq.gz",
        "data/raw_reads/{source}_reads/{id}_2.fastq.gz"
    params:
        source=lambda wc: SAMPLE_SOURCE[wc.id]
    conda:
        "../envs/get_data.yml"
    shell:
        """
        mkdir -p data/raw_reads/{params.source}_reads
        fastq-dump --split-files --gzip \
            -O data/raw_reads/{params.source}_reads \
            {wildcards.id}
        """

