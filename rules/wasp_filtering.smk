rule wasp_filter:
    input:
        bam = "output/aligned/star/{id}.Aligned.sortedByCoord.out.bam",
        snp_dir = "data/genome/known_sites/snp_txt",
        index = "data/index/star"
    output:
        bam = "output/aligned/wasp/{id}.wasp_filtered.bam"
    log:
        "logs/wasp/{id}.log"
    threads: 8
    resources:
        mem_mb = 20000
    conda: "../envs/mapping.yml"
    shell:
        """
        mkdir -p output/aligned/wasp/{wildcards.id}

        python {config[wasp][path]}/find_intersecting_snps.py \
            --is_paired_end \
            --is_sorted \
            --snp_dir {input.snp_dir} \
            --output_dir output/aligned/wasp/{wildcards.id} \
            {input.bam} \
            &> {log}

        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn \
                 output/aligned/wasp/{wildcards.id}/{wildcards.id}.Aligned.sortedByCoord.out.remap.fq1.gz \
                 output/aligned/wasp/{wildcards.id}/{wildcards.id}.Aligned.sortedByCoord.out.remap.fq2.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix output/aligned/wasp/{wildcards.id}/remap_ \
             --outSAMtype BAM SortedByCoordinate \
             &>> {log}

        python {config[wasp][path]}/filter_remapped_reads.py \
            output/aligned/wasp/{wildcards.id}/{wildcards.id}.Aligned.sortedByCoord.out.to.remap.bam \
            output/aligned/wasp/{wildcards.id}/remap_Aligned.sortedByCoord.out.bam \
            output/aligned/wasp/{wildcards.id}/{wildcards.id}.wasp_filtered.unsorted.bam \
            &>> {log}

        # Sort and index the final BAM
        samtools sort -@ {threads} -o {output.bam} \
            output/aligned/wasp/{wildcards.id}/{wildcards.id}.wasp_filtered.unsorted.bam
        samtools index {output.bam}
        """

