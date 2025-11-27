# ============================================================
# Rule: Generate per-sample heterozygous SNP VCF
# ============================================================
rule per_sample_het_vcf:
    input:
        vcf = "output/variants/dna_joint.filtered.het.vcf.gz",
        fasta = "data/genome/hg38.fasta"
    output:
        vcf = "output/variants/per_sample/{sample}.het.vcf.gz"
    log:
        "logs/variants/per_sample_{sample}.log"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk SelectVariants \
            -V {input.vcf} \
            -R {input.fasta} \
            -sn {wildcards.sample} \
            -O {output.vcf}

        # Only index; no extra compression
        tabix -f -p vcf {output.vcf}
        """

