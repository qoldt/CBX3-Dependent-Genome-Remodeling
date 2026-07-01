#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

# Input/output directories
input_dir <- "/Users/kappa/SynologyDrive/MINUTE/peaks/maxgap2000"
output_dir <- "/Users/kappa/SynologyDrive/MINUTE/peaks/maxgap2000/master_peaks_consensus"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Histone marks
marks <- c("H3K9me2", "H3K9me3", "H4K20me3", "H3K36me3")

# Minimum number of replicates required
min_support <- 2

for (mark in marks) {
  
  # Get all files for this mark
  pattern <- paste0(mark, ".*broadPeak$")
  bed_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(bed_files) == 0) {
    message("No files found for ", mark)
    next
  }
  
  message("Processing ", mark, " with ", length(bed_files), " replicates.")
  
  # Read all replicates into GRanges list
  gr_list <- list()
  for (file in bed_files) {
    peakRes <- tryCatch(read.table(file, header = FALSE, fill = TRUE), error = function(e) NULL)
    
    if (is.null(peakRes) || nrow(peakRes) == 0) {
      message("Skipping empty or unreadable file: ", file)
      next
    }
    
    # Check at least 3 columns for BED
    if (ncol(peakRes) < 3) {
      message("Skipping malformed file (less than 3 columns): ", file)
      next
    }
    
    gr <- GRanges(seqnames = peakRes[[1]],
                  ranges = IRanges(start = peakRes[[2]], end = peakRes[[3]]),
                  strand = "*")
    
    # Filter standard chromosomes
    gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
    
    # Skip if nothing remains
    if (length(gr) == 0) {
      message("File has no standard chromosomes after filtering: ", file)
      next
    }
    
    gr_list[[file]] <- gr
  }
  
  # Skip if all files were empty or invalid
  if (length(gr_list) == 0) {
    message("No valid peaks found for ", mark, ", skipping.")
    next
  }
  
  # ✅ Convert list to GRanges by concatenating only valid GRanges
  all_gr <- do.call(c, unname(gr_list))
  
  if (!inherits(all_gr, "GRanges")) {
    stop("Merged object is not GRanges for ", mark)
  }
  
  # Merge all peaks first to get candidate regions
  all_peaks <- GenomicRanges::reduce(all_gr)
  
  # Count support: how many replicates overlap each region
  support_counts <- sapply(gr_list, function(gr) countOverlaps(all_peaks, gr))
  replicate_support <- rowSums(support_counts > 0)
  
  # Filter peaks with enough replicate support
  consensus_peaks <- all_peaks[replicate_support >= min_support]
  
  message("Kept ", length(consensus_peaks), " consensus peaks for ", mark)
  
  # Merge final peaks
  master_peaks <- GenomicRanges::reduce(consensus_peaks)
  
  # Export
  out_file <- file.path(output_dir, paste0(mark, "_consensus_masterPeak.bed"))
  export(master_peaks, out_file, format = "bed")
  
  message("Saved consensus master peaks to: ", out_file)
}

message("✅ Consensus master peak generation completed.")
