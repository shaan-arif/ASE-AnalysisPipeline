rule fastqc:
    input:
        R1 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_1.fastq.gz",
        R2 = lambda wc: f"data/raw_reads/{SAMPLE_SOURCE[wc.id]}_reads/{wc.id}_2.fastq.gz"
    output:
        html = "output/qc/fastqc/{source}/{id}_fastqc.html",
        zip  = "output/qc/fastqc/{source}/{id}_fastqc.zip"
    params:
        source = lambda wc: SAMPLE_SOURCE[wc.id],
        outdir = lambda wc: f"output/qc/fastqc/{SAMPLE_SOURCE[wc.id]}"
    conda:
        "../envs/qc.yml"
    threads: 4
    shell:
        """
        fastqc {input.R1} {input.R2} --outdir {params.outdir} --threads {threads}
        """
  
        
rule multiqc:
    input:
        genomic = expand("output/qc/fastqc/genomic/{id}_fastqc.zip", id=GENOMIC_IDS),
        transcriptomic = expand("output/qc/fastqc/transcriptomic/{id}_fastqc.zip", id=TRANSCRIPTOMIC_IDS)
    output:
        genomic_report = "output/qc/multiqc/genomic/multiqc_report.html",
        genomic_failed = "output/qc/multiqc/genomic/failed_samples.txt",
        transcriptomic_report = "output/qc/multiqc/transcriptomic/multiqc_report.html",
        transcriptomic_failed = "output/qc/multiqc/transcriptomic/failed_samples.txt"
    shell:
        """
        # ▶️ Genomic MultiQC
        mkdir -p output/qc/multiqc/genomic
        multiqc output/qc/fastqc/genomic -o output/qc/multiqc/genomic
        column -t output/qc/multiqc/genomic/multiqc_data/multiqc_fastqc.txt
        awk -F '\\t' 'NR > 1 {{ for(i=1;i<=NF;i++) if($i=="fail" || $i=="warn") {{ print $1; next }} }}' \
            output/qc/multiqc/genomic/multiqc_data/multiqc_fastqc.txt \
            | sort | uniq > {output.genomic_failed}
        cat {output.genomic_failed}

        # ▶️ Transcriptomic MultiQC
        mkdir -p output/qc/multiqc/transcriptomic
        multiqc output/qc/fastqc/transcriptomic -o output/qc/multiqc/transcriptomic
        column -t output/qc/multiqc/transcriptomic/multiqc_data/multiqc_fastqc.txt
        awk -F '\\t' 'NR > 1 {{ for(i=1;i<=NF;i++) if($i=="fail" || $i=="warn") {{ print $1; next }} }}' \
            output/qc/multiqc/transcriptomic/multiqc_data/multiqc_fastqc.txt \
            | sort | uniq > {output.transcriptomic_failed}
        cat {output.transcriptomic_failed}
        """

