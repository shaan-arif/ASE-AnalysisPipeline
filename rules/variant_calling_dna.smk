# ============================================================
# Rule: Call variants per WES sample using GATK HaplotypeCaller
# ============================================================
# Description:
#   This rule runs GATK HaplotypeCaller on each per-sample BQSR-processed WES BAM
#   to generate a genomic VCF (GVCF) containing all sites, including non-variant positions.
#   The GVCF format allows later joint genotyping across multiple samples.
#
# Inputs:
#   - BAM: Deduplicated and BQSR-corrected WES BAM for the sample
#   - Reference FASTA and sequence dictionary
#
# Outputs:
#   - GVCF: Per-sample GVCF file (.g.vcf.gz)
#
# Notes:
#   - Uses `-ERC GVCF` mode for generating GVCFs suitable for joint genotyping.
#   - Logs are written to a sample-specific log file for troubleshooting.
# ============================================================
rule haplotypecaller_dna:
    input:
        bam = "output/bqsr/{id}_dedup_bqsr.bam",
        ref = config["reference"]["fasta"],
        dict = config["reference"]["dict"]  
    output:
        gvcf = "output/variants/gvcf/{id}.g.vcf.gz",
        gvcf_index="output/variants/gvcf/{id}.g.vcf.gz.tbi"
    log:
        "logs/variants/haplotypecaller/{id}.log"
    threads: 8
    resources:
        mem_mb=40000
    conda: "../envs/gatk.yml"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            2> {log}
            
        gatk IndexFeatureFile -I {output.gvcf}
        """

# ============================================================
# Rule: Joint genotyping (combine per-sample GVCFs)
# ============================================================
# Description:
#   Combines all per-sample GVCFs and performs joint genotyping to produce
#   a multi-sample VCF of high-quality variant calls.
#
# Inputs:
#   - List of per-sample GVCFs
#   - Reference FASTA
#
# Outputs:
#   - Joint VCF: Multi-sample VCF with genotypes for all samples
#
# Notes:
#   - First, CombineGVCFs merges individual GVCFs into a single combined GVCF.
#   - Then, GenotypeGVCFs performs variant calling across all samples simultaneously.
#   - This step improves genotyping accuracy, especially for rare variants and low-depth sites.
#   - Log file captures errors or warnings from both steps.
# ============================================================
rule genotype_gvcfs_dna:
    input:
        gvcfs = expand("output/variants/gvcf/{id}.g.vcf.gz", id=GENOMIC_IDS),
        ref = config["reference"]["fasta"]
    output:
        vcf = "output/variants/dna_joint.vcf.gz"
    log:
        "logs/variants/genotype_gvcfs.log"
    threads: 8
    resources:
        mem_mb=40000
    conda: "../envs/gatk.yml"
    params:
        gvcf_list = lambda wildcards, input: " ".join(f"-V {g}" for g in input.gvcfs)
    shell:
        """
        # Step 1: Combine all per-sample GVCFs into a single combined GVCF
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.gvcf_list} \
            -O output/variants/combined.g.vcf.gz

        # Step 2: Joint genotyping across all samples
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V output/variants/combined.g.vcf.gz \
            -O {output.vcf} \
            2> {log}
        """

