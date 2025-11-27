# ============================================================
# Rule: Mark Duplicates
# ============================================================
# This rule marks PCR duplicates in the sorted BAM file.
# It uses GATK's MarkDuplicatesSpark for parallel duplicate marking.
# Input: Sorted BAM from BWA-MEM alignment
# Output: Deduplicated BAM (temporary) + duplication metrics file
# Notes:
# - Temporary BAM is deleted after downstream rules complete.
# - Duplication metrics are stored for QC.
# ============================================================
# ============================================================
# Rule: Mark Duplicates (Non-Spark, Stable)
# ============================================================
# This version uses GATK MarkDuplicates (not Spark).
# It's safer and works well on single-node servers.
# ============================================================
rule mark_duplicates_dna:
    input:
        bam = "output/aligned/bwa/{id}.sorted.bam"
    output:
        bam = "output/aligned/dna_dedup/{id}_dedup.bam",
        bai = "output/aligned/dna_dedup/{id}_dedup.bai",
        metrics = "output/aligned/dna_dedup/{id}_dup_metrics.txt"
    log:
        "logs/markdup_dna/{id}.log"
    threads: 8
    resources:
        mem_mb=16000
    conda: "../envs/gatk.yml"
    shell:
        """
        mkdir -p output/aligned/dna_dedup tmp
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --TMP_DIR tmp \
            2> {log}
        samtools index {output.bam}
        """


# ============================================================
# Rule: Base Recalibration (BQSR)
# ============================================================
# This rule generates a Base Quality Score Recalibration table.
# GATK BaseRecalibrator analyzes patterns of covariation in the sequence dataset.
# Input: Deduplicated BAM + reference genome + known variant sites
# Output: Recalibration table (.table)
# Notes:
# - Known sites (VCFs) are used to model systematic errors.
# - This step is essential to improve variant calling accuracy.
# ============================================================
rule bqsr_dna:
    input:
        bam = "output/aligned/dna_dedup/{id}_dedup.bam",
        ref = config["reference"]["fasta"],
        known_sites = config["dbsnp"]["vcf"]
    output:
        recal_table = "output/bqsr/{id}.recal.table"
    log:
        "logs/bqsr/{id}.log"
    threads: 8
    resources:
        mem_mb=24000
    conda: "../envs/gatk.yml"
    shell:
        """
        gatk BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            --known-sites {input.known_sites} \
            -O {output.recal_table} \
            2> {log}
        """


# ============================================================
# Rule: Apply BQSR
# ============================================================
# This rule applies the recalibration table to the deduplicated BAM.
# Input: Deduplicated BAM + recalibration table + reference genome
# Output: Recalibrated BAM (.bam) + BAM index (.bai)
# Notes:
# - Output BAM is used for downstream variant calling.
# - Indexing is performed immediately after applying BQSR.
# ============================================================
rule apply_bqsr_dna:
    input:
        bam = "output/aligned/dna_dedup/{id}_dedup.bam",
        table = "output/bqsr/{id}.recal.table",
        ref = config["reference"]["fasta"]
    output:
        bam = "output/bqsr/{id}_dedup_bqsr.bam"
    log:
        "logs/apply_bqsr/{id}.log"
    threads: 8
    resources:
        mem_mb=6000
    conda: "../envs/gatk.yml"
    shell:
        """
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} \
            2> {log}
        samtools index {output.bam}
        """

# ============================================================
# Rule: Flagstat after BQSR
# ============================================================
# This rule generates alignment statistics for the final BQSR-applied BAM.
# Input: Recalibrated BAM
# Output: Flagstat summary text file
# Notes:
# - Provides QC metrics like total reads, mapped reads, duplicate reads.
# - Useful for ensuring that BQSR did not negatively affect mapping statistics.
# ============================================================
rule flagstat_dedup_bqsr_dna:
    input:
        bam = "output/bqsr/{id}_dedup_bqsr.bam"
    output:
        stats = "output/bqsr/{id}_dedup_bqsr_flagstat.txt"
    log:
        "logs/bqsr_flagstat/{id}.log"
    threads: 4
    conda: "../envs/mapping.yml"
    shell:
        "samtools flagstat {input.bam} > {output.stats} 2> {log}"

