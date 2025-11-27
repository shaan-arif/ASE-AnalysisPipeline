# ==========================
# Rule: Trimmomatic for transcriptomic data
# ==========================
rule trimmomatic:
    input:
        R1 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_1.fastq.gz",
        R2 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_2.fastq.gz"
    output:
        R1_paired = "data/trimmed_reads/{id}_1.fq.gz",
        R1_unpaired = "data/trimmed_reads/{id}_1U.fq.gz",
        R2_paired = "data/trimmed_reads/{id}_2.fq.gz",
        R2_unpaired = "data/trimmed_reads/{id}_2U.fq.gz",
        log = "data/trimmed_reads/qc/{id}_trimmomatic.log"
    threads: 8
    conda: "../envs/qc.yml"
    params:
        adapters = "tools/trimmomatic/adapters/TruSeq3-PE.fa",
        minlen = 36,
        quality = 20
    shell:
        """
        mkdir -p data/trimmed_reads/qc

        trimmomatic PE -threads {threads} \
            {input.R1} {input.R2} \
            {output.R1_paired} {output.R1_unpaired} \
            {output.R2_paired} {output.R2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:{params.quality} TRAILING:{params.quality} \
            SLIDINGWINDOW:4:{params.quality} MINLEN:{params.minlen} \
            2> {output.log}
        """


# ==========================
# Rule: MultiQC after Trimming
# ==========================
rule multiqc_trimmed:
    input:
        expand("data/trimmed_reads/qc/{id}_trimmomatic.log", id=TRANSCRIPTOMIC_IDS)
    output:
        html = "data/trimmed_reads/qc/multiqc_report.html"
    threads: 2
    conda: "../envs/qc.yml"
    shell:
        """
        mkdir -p data/trimmed_reads/qc

        multiqc data/trimmed_reads/qc/ \
            --filename {output.html} \
            --title "Trimmomatic QC Summary" \
            --quiet
        """

