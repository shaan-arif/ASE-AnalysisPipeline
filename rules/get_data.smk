# ==========================
# Rule: Download SRA Runs from BioProject PRJNA675090
# ==========================
rule get_data:
    output:
        "data/raw_reads/{source}_reads/{id}_1.fastq.gz",
        "data/raw_reads/{source}_reads/{id}_2.fastq.gz"
    params:
        source=lambda wc: SAMPLES.loc[SAMPLES["Run"] == wc.id, "LibrarySource"].values[0].lower()
    conda:
        "../envs/get_data.yml"
    shell:
        """
        fastq-dump --split-files --gzip \
            -O data/raw_reads/{params.source}_reads \
            {wildcards.id}
        """

