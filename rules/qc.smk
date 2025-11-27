# ============================================================
# Rule: fastqc
# Purpose: Perform quality control on raw paired-end reads using FastQC.
# Generates both HTML reports and zipped output files for each sample.
# ============================================================
rule fastqc:
    input:
        R1 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_1.fastq.gz",
        R2 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_2.fastq.gz"
    output:
        html1 = "output/qc/fastqc/{source}/{id}_1_fastqc.html",
        zip1  = "output/qc/fastqc/{source}/{id}_1_fastqc.zip",
        html2 = "output/qc/fastqc/{source}/{id}_2_fastqc.html",
        zip2  = "output/qc/fastqc/{source}/{id}_2_fastqc.zip"
    params:
        source = lambda wc: SAMPLE_SOURCE[wc.id],
        outdir = lambda wc: f"output/qc/fastqc/{SAMPLE_SOURCE[wc.id]}"
    threads: 4
    conda: "../envs/qc.yml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.R1} {input.R2} --outdir {params.outdir} --threads {threads}
        """



# ============================================================
# Rule: multiqc
# Purpose: Aggregate FastQC results across all samples into a single
# interactive MultiQC report per sample type (genomic and transcriptomic).
# Also generates a list of failed/warning samples for quick inspection.
# ============================================================
rule multiqc:
    input:       
        genomic = expand("output/qc/fastqc/genomic/{id}_1_fastqc.zip", id=GENOMIC_IDS),
        transcriptomic = expand("output/qc/fastqc/transcriptomic/{id}_1_fastqc.zip", id=TRANSCRIPTOMIC_IDS)
        
    output:        
        genomic_report = "output/qc/multiqc/genomic/multiqc_report.html",
        genomic_failed = "output/qc/multiqc/genomic/failed_samples.txt",
        transcriptomic_report = "output/qc/multiqc/transcriptomic/multiqc_report.html",
        transcriptomic_failed = "output/qc/multiqc/transcriptomic/failed_samples.txt"
   
    shell:
        """
        mkdir -p output/qc/multiqc/genomic
        multiqc output/qc/fastqc/genomic -o output/qc/multiqc/genomic

        awk -F '\\t' 'NR>1{{for(i=1;i<=NF;i++) if($i=="fail"||$i=="warn"){{print $1;next}}}}' \
            output/qc/multiqc/genomic/multiqc_data/multiqc_fastqc.txt | sort | uniq > {output.genomic_failed}

        mkdir -p output/qc/multiqc/transcriptomic
        multiqc output/qc/fastqc/transcriptomic -o output/qc/multiqc/transcriptomic
        awk -F '\\t' 'NR>1{{for(i=1;i<=NF;i++) if($i=="fail"||$i=="warn"){{print $1;next}}}}' \
            output/qc/multiqc/transcriptomic/multiqc_data/multiqc_fastqc.txt | sort | uniq > {output.transcriptomic_failed}
        """

