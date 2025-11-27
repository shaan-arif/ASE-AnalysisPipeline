
#===========================================================
#              ASE Read Counting       
#===========================================================

rule ase_counting:
    input:
        bam=lambda wildcards: f"output/aligned/star/{patient_dict[wildcards.patient][0]}.Aligned.sortedByCoord.out.bam",
        vcf=lambda wildcards: f"output/variants/per_sample/{patient_dict[wildcards.patient][1]}.het.vcf.gz",
        ref=config['reference']['fasta']
    output:
        table="output/ase/{patient}.ase_counts.table"
    threads: 4
    conda: "../envs/gatk.yml"
    shell:
        """
        mkdir -p output/ase
        gatk ASEReadCounter \
            -R {input.ref} \
            -I {input.bam} \
            -V {input.vcf} \
            -O {output.table} \
            --min-depth 10 --min-mapping-quality 10 --min-base-quality 20
        """

