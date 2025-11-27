# ============================================================
# Rule: Filter DNA SNPs for ASE analysis with per-sample stats
# ============================================================

rule filter_dna_for_ase:
    input:
        vcf = "output/variants/dna_joint.vcf.gz",
        fasta = "data/genome/hg38.fasta",
        repeatmasker = "others/repeatmasker/hg38_rmsk.sorted.bed"
    output:
        vcf = "output/variants/dna_joint.filtered.het.vcf.gz",
        stats = "output/variants/ase_snp_stats.tsv",
        per_sample_stats = "output/variants/ase_snp_per_sample.tsv"
    log:
        "logs/variants/filter_dna_for_ase.log"
    threads: 8
    resources:
        mem_mb = 32000
    conda: "../envs/gatk.yml"
    shell:
        r"""
        mkdir -p output/variants/temp_files

        # Step 1: Select SNPs only
        echo "[INFO] Step 1: Selecting SNPs only..." | tee -a {log}
        gatk SelectVariants -R {input.fasta} -V {input.vcf} --select-type-to-include SNP \
            -O output/variants/temp_files/temp_snp.vcf.gz

        # Step 2: Keep biallelic SNPs
        echo "[INFO] Step 2: Keeping only biallelic SNPs..." | tee -a {log}
        gatk SelectVariants -R {input.fasta} -V output/variants/temp_files/temp_snp.vcf.gz \
            --restrict-alleles-to BIALLELIC \
            -O output/variants/temp_files/temp_snp_biallelic.vcf.gz

        # Step 3: Keep heterozygous genotypes (DP>=10, GQ>=20)
        echo "[INFO] Step 3: Keeping heterozygous genotypes (DP>=10, GQ>=20)..." | tee -a {log}
        bcftools view -g het -i 'FMT/DP>=10 & FMT/GQ>=20' \
            output/variants/temp_files/temp_snp_biallelic.vcf.gz \
            -Oz -o output/variants/temp_files/temp_het.vcf.gz
        tabix -p vcf output/variants/temp_files/temp_het.vcf.gz

        # Step 3b: Count heterozygous SNPs before filtering
        num_het_before=$(bcftools view -H output/variants/temp_files/temp_het.vcf.gz | wc -l)
        echo "[INFO] Number of heterozygous SNPs before filtering: $num_het_before" | tee -a {log}

        # Step 4: Filter low-quality variants
        echo "[INFO] Step 4: Filtering low-quality variants..." | tee -a {log}
        gatk VariantFiltration -R {input.fasta} -V output/variants/temp_files/temp_het.vcf.gz \
            --filter-expression "QD < 2.0 || FS > 30.0 || MQ < 40.0 || SOR > 3.0" \
            --filter-name "low_quality_variant" \
            -O output/variants/temp_files/temp_filtered_het.vcf.gz

        # Step 5: Keep only PASS variants
        echo "[INFO] Step 5: Keeping only PASS variants..." | tee -a {log}
        bcftools view -f PASS output/variants/temp_files/temp_filtered_het.vcf.gz \
            -Oz -o output/variants/temp_files/temp_het.pass.vcf.gz
        tabix -p vcf output/variants/temp_files/temp_het.pass.vcf.gz

        # Step 6: Remove RepeatMasker regions
        echo "[INFO] Step 6: Removing RepeatMasker regions..." | tee -a {log}
        bcftools view -T ^{input.repeatmasker} \
            output/variants/temp_files/temp_het.pass.vcf.gz -Oz -o output/variants/temp_files/temp_final.vcf.gz
        tabix -p vcf output/variants/temp_files/temp_final.vcf.gz

        # Step 7: Count heterozygous SNPs after all filters
        num_het_after=$(bcftools view -H output/variants/temp_files/temp_final.vcf.gz | wc -l)
        echo "[INFO] Number of heterozygous SNPs after filtering: $num_het_after" | tee -a {log}

        # Step 8: Compress and index final VCF
        cp output/variants/temp_files/temp_final.vcf.gz {output.vcf}
        tabix -p vcf {output.vcf}

        # Step 9: Save joint stats
        echo -e "Stage\tNum_Het_SNPs" > {output.stats}
        echo -e "Before_filters\t$num_het_before" >> {output.stats}
        echo -e "After_filters\t$num_het_after" >> {output.stats}

        # Step 10: Per-sample heterozygous SNP counts
        echo "[INFO] Counting heterozygous SNPs per sample..." | tee -a {log}
        echo -e "Sample\tNum_Het_SNPs_Before\tNum_Het_SNPs_After" > {output.per_sample_stats}

        for sample in $(bcftools query -l output/variants/temp_files/temp_het.vcf.gz); do
            het_before=$(bcftools view -H -s $sample output/variants/temp_files/temp_het.vcf.gz | wc -l)
            het_after=$(bcftools view -H -s $sample output/variants/temp_files/temp_final.vcf.gz | wc -l)
            echo -e "$sample\t$het_before\t$het_after" >> {output.per_sample_stats}
        done

        echo "[SUCCESS] ASE-ready VCF and per-sample stats generated." | tee -a {log}
        """

