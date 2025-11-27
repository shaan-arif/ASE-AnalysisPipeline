# ==========================
# Rule: STAR genome index
# ==========================
# Creates STAR index files from the reference genome FASTA and GTF annotation.
# These indexes are stored in "data/index/star/" and are required before alignment.
# The sjdbOverhang parameter should be set to (read_length - 1).
# Here it is set to 100, suitable for ~101 bp Illumina reads.

rule star_index:
    input:
        fasta = config["reference"]["fasta"],
        gtf   = config["annotation"]["gtf"]
    output:
        directory("data/index/star")
    threads: 8
    resources:
        mem_mb=25600
    conda: "../envs/mapping.yml"
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100
        """


# ==========================
# Rule: STAR alignment
# ==========================
# Aligns paired-end RNA-seq reads (trimmed) to the reference genome using STAR.
# Uses the pre-built STAR index from "data/index/star".
# Produces a coordinate-sorted BAM file, suitable for downstream analyses 
# such as transcript quantification or variant calling.
rule star_align:
    input:
        r1 = "data/trimmed_reads/{id}_1.fq.gz",
        r2 = "data/trimmed_reads/{id}_2.fq.gz",
        index = "data/index/star"
    output:
        bam = "output/aligned/star/{id}.Aligned.sortedByCoord.out.bam",
        stats = "output/aligned/star/{id}_flagstat.txt"
    log:
        "logs/star/{id}.log"
    threads: 8
    resources:
        mem_mb=20000
    conda: "../envs/mapping.yml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrRGline ID:{wildcards.id} SM:{wildcards.id} PL:ILLUMINA LB:lib1 \
             --outFilterMultimapNmax 1 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000 \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --quantMode GeneCounts \
             --outFileNamePrefix output/aligned/star/{wildcards.id}.

        
        # Generate alignment statistics
        samtools flagstat {output.bam} > {output.stats}
        """

