#!/usr/bin/env Rscript

library(data.table)
library(MBASED)
library(GenomicRanges)
library(SummarizedExperiment)
library(BiocParallel)

### ================================
### User paths
### ================================
ase_dir    <- "output/ase"
map_dir    <- "data/annotation/gene_mapping"
mapping_csv <- "config/patient_sample_map.csv"
outdir     <- file.path(ase_dir, "MBASED_results")

min_total_depth <- 10
min_minor_count <- 3
numSim <- 1e5
BPPARAM <- SerialParam()   # change to MulticoreParam(workers=8) if needed

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "rds"), showWarnings = FALSE)
dir.create(file.path(outdir, "summary_tables"), showWarnings = FALSE)
dir.create(file.path(outdir, "input_tables"), showWarnings = FALSE)

### ================================
### Load patient → WES_ID mapping
### ================================
map_samples <- fread(mapping_csv)
setnames(map_samples, c("Patient","RNAseq_ID","WES_ID"))

### ================================
### List ASE files
### ================================
ase_files <- list.files(ase_dir, pattern = "HCM.*ase_counts.table", full.names = TRUE)
if (length(ase_files) == 0) stop("No ASE tables found!")

### ================================
### Main loop
### ================================
all_gene_results <- list()

for (file in ase_files) {
  
  sample_name <- gsub(".ase_counts.table", "", basename(file))
  
  message("\n========== Processing sample: ", sample_name, " ==========")
  
  ### Find matching WES_ID
  wes_id <- map_samples[Patient == sample_name, WES_ID]
  if (length(wes_id) == 0) {
    message("  No WES ID for sample ", sample_name, " — skipping")
    next
  }
  
  map_file <- file.path(map_dir, paste0(wes_id, ".snp_gene_map.tsv"))
  if (!file.exists(map_file)) {
    message("  Gene map missing: ", map_file, " — skipping")
    next
  }
  
  ### Load ASE
  ase <- fread(file)
  ase[, position := as.numeric(position)]
  
  ### Load SNP → gene map
  gmap <- fread(map_file, header = FALSE,
                col.names = c("contig", "position", "gene"))
  gmap[, position := as.numeric(position)]
  
  ### Merge
  dt <- merge(ase, gmap, by = c("contig","position"), all.x = TRUE)
  dt <- dt[!is.na(gene)]
  
  ### Apply filters
  dt <- dt[totalCount >= min_total_depth &
             pmin(refCount, altCount) >= min_minor_count]
  
  if (nrow(dt) == 0) {
    message("  No SNPs passing filters — skipping")
    next
  }
  
  ### Build MBASED input
  mdf <- data.frame(
    chr = dt$contig,
    pos = as.integer(dt$position),
    refAllele = dt$refAllele,
    altAllele = dt$altAllele,
    refCount = dt$refCount,
    altCount = dt$altCount,
    gene = dt$gene,
    stringsAsFactors = FALSE
  )
  
  fwrite(as.data.table(mdf),
         file.path(outdir, "input_tables",
                   paste0(sample_name, "_MBASED_input.tsv")),
         sep = "\t")
  
  lociAllele1Counts <- matrix(mdf$refCount, ncol = 1,
                              dimnames = list(NULL, sample_name))
  lociAllele2Counts <- matrix(mdf$altCount, ncol = 1,
                              dimnames = list(NULL, sample_name))
  lociAllele1CountsNoASEProbs <- matrix(0.5, nrow = nrow(mdf), ncol = 1)
  lociCountsDispersions <- matrix(0, nrow = nrow(mdf), ncol = 1)
  
  rr <- GRanges(seqnames = mdf$chr,
                ranges = IRanges(start = mdf$pos, end = mdf$pos),
                aseID = mdf$gene)
  
  ASESummarizedExperiment <- SummarizedExperiment(
    assays = list(
      lociAllele1Counts = lociAllele1Counts,
      lociAllele2Counts = lociAllele2Counts,
      lociAllele1CountsNoASEProbs = lociAllele1CountsNoASEProbs,
      lociCountsDispersions = lociCountsDispersions
    ),
    rowRanges = rr
  )
  
  ### Run MBASED
  res <- tryCatch({
    runMBASED(
      ASESummarizedExperiment = ASESummarizedExperiment,
      isPhased = FALSE,
      numSim = numSim,
      BPPARAM = BPPARAM
    )
  }, error = function(e) {
    message(" MBASED failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(res)) next
  
  saveRDS(res,
          file = file.path(outdir, "rds",
                           paste0(sample_name, "_MBASED.rds")))
  
  ### Extract gene-level results
  maf <- assays(res)[["majorAlleleFrequency"]][,1]
  pval <- assays(res)[["pValueASE"]][,1]
  het  <- assays(res)[["pValueHeterogeneity"]][,1]
  genes <- rownames(assays(res)[[1]])
  
  out <- data.table(
    gene = genes,
    majorAlleleProportion = maf,
    pValue = pval,
    pValueHeterogeneity = het,
    sample = sample_name
  )
  
  fwrite(out,
         file.path(outdir, "summary_tables",
                   paste0(sample_name, "_MBASED_gene_summary.tsv")),
         sep = "\t")
  
  all_gene_results[[sample_name]] <- out
  
  message("  Done sample ", sample_name)
}

### ================================
### Combine all samples
### ================================
combined <- rbindlist(all_gene_results, fill = TRUE)
fwrite(combined, file.path(outdir, "gene_level_ASE_MBASED_all_samples.tsv"), sep = "\t")

message("\nAll MBASED processing completed.")
