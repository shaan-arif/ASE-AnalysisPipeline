
configfile: "config/config.yaml"

import pandas as pd

chroms = [str(i) for i in range(1,23)] + ["X","Y","MT"]
snp_txts = expand("data/genome/known_sites/snp_txts/{chr}.snps.txt.gz", chr=chroms)


# Load mapping table
sample_map = pd.read_csv("config/patient_sample_map.csv", sep="\t")

# Create a dictionary: Patient -> (RNAseq_ID, WES_ID)
patient_dict = {row['Patient']:(row['RNAseq_ID'], row['WES_ID']) for _, row in sample_map.iterrows()}



# Load sample table and extract Run IDs
SAMPLES = pd.read_csv(config["sample_table"])
SRA_IDS = SAMPLES["Run"].tolist()
SAMPLE_SOURCE = dict(zip(SAMPLES["Run"], SAMPLES["LibrarySource"].str.lower()))
GENOMIC_IDS = SAMPLES.loc[SAMPLES["LibrarySource"].str.lower() == "genomic", "Run"].tolist()
TRANSCRIPTOMIC_IDS = SAMPLES.loc[SAMPLES["LibrarySource"].str.lower() == "transcriptomic", "Run"].tolist()

include: "rules/get_reference.smk"
include: "rules/get_data.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"         
include: "rules/wasp_filtering.smk"
include: "rules/dna_mapping.smk" 
include: "rules/dedup_bqsr.smk"
include: "rules/variant_calling_dna.smk"  
include: "rules/variant_filtering_dna.smk"    
include: "rules/per_sample_het_VCF.smk"                   
include: "rules/multiqc_final.smk"
include: "rules/ASEReadcounter.smk" 
include: "rules/MBASED.smk"
include:"rules/plots.smk"



rule all:
    input:
        #----------------------------
        #   MultiQC report
        #----------------------------
        "output/qc/multiqc/genomic/multiqc_report.html",
        "output/qc/multiqc/genomic/failed_samples.txt",
        "output/qc/multiqc/transcriptomic/multiqc_report.html",
        "output/qc/multiqc/transcriptomic/failed_samples.txt",
        # --------------------------
        # Trimmed transcriptomic reads
        # --------------------------
        expand("data/trimmed_reads/{id}_1.fq.gz", id=TRANSCRIPTOMIC_IDS),
        expand("data/trimmed_reads/{id}_2.fq.gz", id=TRANSCRIPTOMIC_IDS),
        # --------------------------
        # Genomic mapping (BWA → sorted BAM)
        # --------------------------
        expand("output/aligned/bwa/{id}_flagstat.txt", id=GENOMIC_IDS),
        expand("output/aligned/dna_dedup/{id}_dup_metrics.txt", id=GENOMIC_IDS),
        expand("output/bqsr/{id}_dedup_bqsr_flagstat.txt", id=GENOMIC_IDS),
        # --------------------------
        # Transcriptomic mapping (STAR → sorted BAM)
        # --------------------------
        expand("output/aligned/star/{id}_flagstat.txt", id=TRANSCRIPTOMIC_IDS),
        expand("output/aligned/wasp/{id}.wasp_filtered.bam", id=TRANSCRIPTOMIC_IDS),      
        # --------------------------
        # Final MultiQC across DNA + RNA
        # --------------------------
        "output/qc/multiqc/final/multiqc_report.html",
        # ------------------------------
        #        Final Filtered joint vcf
        #-------------------------------
        #"output/variants/dna_joint.filtered.het.vcf.gz",
        #--------------------------------------
        #        Per Sample Het VCF
        #--------------------------------------
        expand("output/variants/per_sample/{id}.het.vcf.gz", id=GENOMIC_IDS),      
        #--------------------------------------
        #     ASE Read Counts Patient Wise
        #--------------------------------------
        expand("output/ase/{patient}.ase_counts.table", patient=list(patient_dict.keys())),
        "output/ase/MBASED_results/plots/done.flag",
