# -----------------------------------------------------------
#                  MBASED Analysis
# -----------------------------------------------------------

rule m_based:
    input:
        ase_tables = expand("output/ase/{patient}.ase_counts.table",
                            patient=list(patient_dict.keys())),
        map_csv = "config/patient_sample_map.csv",
        gene_map_dir = "data/annotation/gene_mapping",
        script = "scripts/MBASED.R"
    output:
        combined = "output/ase/MBASED_results/gene_level_ASE_MBASED_all_samples.tsv"
    conda:
        "envs/r_mbased.yml"
    shell:
        """
        Rscript {input.script}
        """

