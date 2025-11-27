# ============================================================
# Rule: Download Reference Genome
# ============================================================
rule download_reference:
    output:
        fasta = config["reference"]["fasta"],
        fai   = config["reference"]["fasta"] + ".fai",
        dict  = config["reference"]["dict"]
    params:
        url = config["reference"]["ftp"]
    log:
        "logs/get_data/download_reference.log"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta}.gz {params.url} 2> {log}
        gunzip -c {output.fasta}.gz > {output.fasta}
        rm {output.fasta}.gz
        
         # Create FASTA index (.fai)
        samtools faidx {output.fasta}
        
        gatk CreateSequenceDictionary \
        -R {output.fasta} \
        -O {output.dict}

        """

# ============================================================
# Rule: Download Gene Annotation (GTF)
# ============================================================
rule download_annotation:
    output:
        gtf = config["annotation"]["gtf"]
    params:
        url = config["annotation"]["ftp"]
    log:
        "logs/get_data/download_annotation.log"
    shell:
        """
        mkdir -p $(dirname {output.gtf})
        wget -O {output.gtf}.gz {params.url} 2> {log}
        gunzip -c {output.gtf}.gz > {output.gtf}
        rm {output.gtf}.gz
        """

# ============================================================
# Rule: Download dbSNP (VCF + Index)
# ============================================================

#rule download_ensembl_vcf:
#    output:
#        vcf = "data/genome/ensembl_vcf/homo_sapiens-chr{chr}.vcf.gz",
#        index = "data/genome/ensembl_vcf/homo_sapiens-chr{chr}.vcf.gz.tbi"
#    params:
#        ftp = "https://ftp.ensembl.org/pub/release-111/variation/vcf/homo_sapiens/homo_sapiens-chr{chr}.vcf.gz"
#    shell:
#        """
#        mkdir -p data/genome/ensembl_vcf 
#        wget -O {output.vcf} {params.ftp}
#        tabix -p vcf {output.vcf}
#        """
#
# Then you can expand it:
#vcfs = expand("data/genome/ensembl_vcf/homo_sapiens-chr{chr}.vcf.gz", chr=chroms)

# ============================================================
# Rule: Download dbSNP VCF + Index
# ============================================================
rule download_dbsnp:
    output:
        vcf = config["dbsnp"]["vcf"],
        index = config["dbsnp"]["index"]
    params:
        url = config["dbsnp"]["ftp"]
    log:
        "logs/get_data/download_dbsnp.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        if [ ! -f {output.vcf} ]; then
            wget -O {output.vcf} {params.url} 2> {log}
            tabix -p vcf {output.vcf}
        fi
        """

# ============================================================
# Rule: Download RepeatMasker Annotation
# ============================================================
rule download_repeatmasker:
    output:
        rm = config["repeatmasker"]["rm"]
    params:
        url = config["repeatmasker"]["ftp"]
    log:
        "logs/get_data/download_repeatmasker.log"
    shell:
        """
        mkdir -p $(dirname {output.rm})
        wget -O {output.rm} {params.url} 2> {log}       
        """


# ============================================================
# Rule: Convert RepeatMasker .out.gz to sorted BED
# ============================================================
# ============================================================
# Rule: Convert RepeatMasker .out.gz to sorted BED
# ============================================================
rule convert_repeatmasker_to_bed:
    input:
        rm = config["repeatmasker"]["rm"]
    output:
        bed = "others/repeatmasker/hg38_rmsk.sorted.bed"
    log:
        "logs/get_data/convert_repeatmasker_to_bed.log"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        gunzip -c {input.rm} | \
        grep -v "^   SW" | \
        grep -v "^score" | \
        awk 'NR>3 {{print $5 "\\t" $6-1 "\\t" $7 "\\t" $10}}' | \
        sort -k1,1 -k2,2n > {output.bed}
        """

