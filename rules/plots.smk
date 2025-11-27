# -----------------------------------------------------------
#                 Generate ASE Plots
# -----------------------------------------------------------

rule ase_plots:
    input:
        combined = "output/ase/MBASED_results/gene_level_ASE_MBASED_all_samples.tsv",
        script = "scripts/plots.R"
    output:
        flag = "output/ase/MBASED_results/plots/done.flag"
    
    shell:
        """
        mkdir -p output/ase/MBASED_results/plots
        Rscript {input.script}
        touch {output.flag}
        """

