# ======================================
# Rule: Index reference genome for BWA
# Generates: .amb, .ann, .bwt, .pac, .sa index files.
# =======================================
rule bwa_index:
    input:
        ref = config['reference']['fasta']
    output:
        expand("data/index/bwa/hg38.{ext}",
               ext=["amb", "ann", "bwt", "pac", "sa"])
    conda:
        "../envs/mapping.yml"
    threads: 4
    resources:
        mem_mb=25600
    shell:
        """
        # Indexing reference genome
        bwa index -p data/index/bwa/hg38 {input.ref}
        """


# ===================================
# Rule: BWA-MEM + sort + flagstat
# Purpose: Align paired-end genomic reads to reference genome using BWA-MEM,
# sort the resulting BAM file, and compute alignment statistics.
# Generates: Sorted BAM and flagstat text file per sample.
# ===================================
rule bwa_mem:
    input:
        ref = config['reference']['fasta'],
        r1 = "data/raw_reads/genomic_reads/{id}_1.fastq.gz",
        r2 = "data/raw_reads/genomic_reads/{id}_2.fastq.gz",
        index = expand("data/index/bwa/hg38.{ext}",
                       ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        bam = "output/aligned/bwa/{id}.sorted.bam",
        stats = "output/aligned/bwa/{id}_flagstat.txt"
    log:
        "logs/bwa_mem/{id}.log"
    threads: 8
    resources:
        mem_mb=20000
    conda: "../envs/mapping.yml"
    shell:
        r"""
        set -euo pipefail

        # BWA-MEM alignment + sort in one go
        bwa mem -t {threads} \
            -R '@RG\tID:{wildcards.id}\tSM:{wildcards.id}\tPL:ILLUMINA\tLB:lib1' \
            data/index/bwa/hg38 {input.r1} {input.r2} 2> {log} | \
        samtools view -Sb - | \
        samtools sort -@ {threads} -o {output.bam} -

        # Alignment QC
        samtools flagstat {output.bam} > {output.stats}
        """
