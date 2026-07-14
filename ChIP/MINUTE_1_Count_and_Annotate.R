###
###
### This script takes scaled bigwigs from MINUTE Multiplex ChIPseq and called peaks as inputs
### Given the peak files it queries the ChIP signal in the peak file range and then returns counts
### Signal is tested between genotypes using DESeq2, without the scaling step, as these bigwigs are scaled already
### All regions are intersected with Genomic annotation, gene annotation and repeat annotation.
### It then returns a list of annotated dataframes for each ChIPseq

source("config.R")

# Generate and store count matrices

count_tables <- list()  # in-memory storage

for (mark in marks) {
  cat("\n=== Generating counts for:", mark, "===\n")

  region_file <- regions[[mark]]
  bigwigs <- get_bigwig_files(mark)
  
  if (!file.exists(region_file)) {
    cat("⚠️  Skipping", mark, "- region file not found:", region_file, "\n")
    next
  }
  
  peaks <- import(region_file)
  seqlevelsStyle(peaks) <- "NCBI"
  cat(" - Imported", length(peaks), "regions\n")
  
  results <- mclapply(seq_along(bigwigs), function(i) {
    bw_file <- bigwigs[[i]]
    sample_name <- names(bigwigs)[i]
    
    if (!file.exists(bw_file)) {
      message(sprintf("Missing: %s", bw_file))
      return(rep(NA_real_, length(peaks)))
    }
    
    message(sprintf("[%d/%d] %s", i, length(bigwigs), sample_name))
    tryCatch({
      vals <- bw_loci(bw_file, peaks)
      mcols(vals)[[1]]
    }, error = function(e) {
      message(sprintf("Error in %s: %s", sample_name, e$message))
      rep(NA_real_, length(peaks))
    })
  }, mc.cores = detectCores(), mc.preschedule = FALSE)
  
  # Clean and build matrix
  results_clean <- lapply(results, function(x) {
    if (is.numeric(x) && length(x) == length(peaks)) x else rep(NA_real_, length(peaks))
  })
  count_matrix <- do.call(cbind, results_clean)
  colnames(count_matrix) <- names(bigwigs)
  row_ids <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
  rownames(count_matrix) <- row_ids
  count_df <- as.data.frame(count_matrix)
  
  # Save to memory and file
  count_tables[[mark]] <- count_df
  out_file <- file.path(counts_dir, paste0(mark, "_bigwig_counts.tsv"))
  fwrite(cbind(peak_id = row_ids, count_df), out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("✅ Saved:", out_file, "\n")
}

###
# DESeq2 Loop over each ChIP
###

dds_objects <- list()
norm_counts <- list()
deseq_results <- list()

for (mark in names(count_tables)) {
  cat("\n=== Running DESeq2 for:", mark, "===\n")
  
  raw_df <- count_tables[[mark]]
  
  # Rebuild rownames if chr/start/end present
  if (all(c("chr", "start", "end") %in% names(raw_df))) {
    coord_strings <- paste(raw_df$chr, raw_df$start, raw_df$end, sep = ":")
    raw_counts <- raw_df[, !(names(raw_df) %in% c("chr", "start", "end"))]
    stopifnot(length(coord_strings) == nrow(raw_counts))
    rownames(raw_counts) <- coord_strings
  } else {
    raw_counts <- raw_df
    coord_strings <- rownames(raw_df)
  }
  
  raw_counts_int <- round(raw_counts)  # Ensure integer counts for DESeq2
  
  # Filter low-count rows
  keep <- which(rowSums(raw_counts_int) > 10)
  filtered_counts <- raw_counts_int[keep, ]
  
  # Sample metadata — genotype pulled from the sample sheet by column name,
  # so it never drifts from the actual bigWig set. WT is the reference level.
  geno <- samples$genotype[match(colnames(filtered_counts), samples$sample_id)]
  metadata <- data.frame(
    row.names = colnames(filtered_counts),
    condition = factor(as.character(geno), levels = c("WT", "HP1gKO"))
  )
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~ condition)
  sizeFactors(dds) <- rep(1, ncol(filtered_counts))  # Use scaled bigWig input SO NO SCALING
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  # Save outputs
  dds_objects[[mark]] <- dds
  norm_counts[[mark]] <- counts(dds, normalized = TRUE)
  
  res <- as.data.frame(results(dds, contrast = c("condition", "HP1gKO", "WT")))
  rownames(res) <- rownames(filtered_counts)
  deseq_results[[mark]] <- res
  
  cat("✅ DESeq2 complete for:", mark, "\n")
}

# Plot distributions of Log2FoldChange by ChIP
par(mfrow = c(2, 3))  # Adjust as needed

for (chip in names(deseq_results)) {
  res <- deseq_results[[chip]]
  res <- res[!is.na(res$log2FoldChange), ]
  
  hist(res$log2FoldChange,
       breaks = 100,
       main = chip,
       xlab = "log2FC",
       col = "gray80",
       border = "white")
}


# Plot distributions of pvalues by ChIP
par(mfrow = c(2, 3))  # Adjust as needed

for (chip in names(deseq_results)) {
  res <- deseq_results[[chip]]
  res <- res[!is.na(res$pvalue), ]
  
  hist(res$pvalue,
       breaks = 100,
       main = chip,
       xlab = "pvalue",
       col = "gray80",
       border = "white")
}



### LOAD ANNOTATIONS

# --- Load Repeat + TAD annotation (shared loader from config.R) ---
ann     <- load_annotation()
line_gr <- ann$line
sine_gr <- ann$sine
ltr_gr  <- ann$ltr
tad_gr  <- ann$tad


# --- Annotation Loop ---

annotated_results <- list()

for (mark in names(deseq_results)) {
  cat("\n=== Annotating:", mark, "===\n")
  
  res <- deseq_results[[mark]]
  res <- res[complete.cases(res), ]  # Remove NAs from DESeq
  
  # Extract coordinates from rownames
  peak_coords <- rownames(res)
  coord_parts <- do.call(rbind, strsplit(peak_coords, ":|-"))
  colnames(coord_parts) <- c("chr", "start", "end")
  
  res$chr <- coord_parts[, 1]
  res$start <- as.numeric(coord_parts[, 2])
  res$end <- as.numeric(coord_parts[, 3])
  
  peak_gr <- GRanges(seqnames = res$chr, ranges = IRanges(res$start, res$end))
  seqlevelsStyle(peak_gr) <- "UCSC"
  
  # Genomic annotation
  # TO DO : assess priority, very large domains get labelled as Promoter. 
  # Perhaps a custom annotation for large domains can be added downstream
  peak_anno <- annotatePeak(peak_gr,
                            TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                            tssRegion = c(-1500, 1500),
                            verbose = FALSE)
  anno_df <- as.data.frame(peak_anno)
  res$genomic_region <- anno_df$annotation
  res$genomic_region <- sub("\\s*\\(.*\\)", "", res$genomic_region)
  res$gene_id <- anno_df$geneId
  # This is also a point of concern, some Distal Intergenic Regions have extremely large distance to TSS, 
  # So one should be careful when annotating with gene name as several regions are very far from the gene.
  res$distance_to_tss <- anno_df$distanceToTSS
  
  # TAD boundaries
  tad_anno <- subsetByOverlaps(peak_gr, tad_gr)
  seqlevelsStyle(tad_anno) <- "NCBI" # because chr on res is in this style
  tad_df <- as.data.frame(tad_anno)
  colnames(tad_df)[1] <- "chr"
  tad_df <- tad_df[1:3] # keep only coordinate columns (chr, start, end)
  tad_df$overlaps_with_Tad_boundary <- TRUE
  
  res <- left_join(res, tad_df, by = c("chr","start","end"))
  res$overlaps_with_Tad_boundary[is.na(res$overlaps_with_Tad_boundary)] <- FALSE
  
  # Gene SYMBOLs
  gene_annot <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(res$gene_id),
    columns = c("SYMBOL", "ENSEMBL"),
    keytype = "ENTREZID"
  )
  gene_annot <- gene_annot[!duplicated(gene_annot$ENTREZID), ]
  res <- left_join(res, gene_annot, by = c("gene_id" = "ENTREZID"))
  
  # Label: gene or coordinate (It seems that currently most things get a gene label)
  res$label <- ifelse(!is.na(res$SYMBOL),
                      res$SYMBOL,
                      paste0(res$chr, ":", res$start, "-", res$end))
  
  # Repeat classification
  
  
  # Ensure seqlevel styles match
  seqlevelsStyle(line_gr) <- seqlevelsStyle(peak_gr)
  seqlevelsStyle(sine_gr) <- seqlevelsStyle(peak_gr)
  seqlevelsStyle(ltr_gr)  <- seqlevelsStyle(peak_gr)
  
  # Combine all repeats into one GRanges and keep metadata
  repeat_gr <- c(line_gr, sine_gr, ltr_gr)
  
  # Find overlaps once
  hits <- findOverlaps(peak_gr, repeat_gr, ignore.strand = TRUE)
  
  # Default values (no overlaps)
  res$repeat_class  <- rep("none", length(peak_gr))
  res$repeat_name   <- rep(NA_character_, length(peak_gr))
  res$repeat_family <- rep(NA_character_, length(peak_gr))
  
  if (length(hits) > 0) {
    df_hits <- data.table(
      q         = queryHits(hits),                              # index of peak_gr
      repClass  = mcols(repeat_gr)$repClass[subjectHits(hits)],
      repName   = mcols(repeat_gr)$repName[subjectHits(hits)],
      repFamily = mcols(repeat_gr)$repFamily[subjectHits(hits)]
    )
    
    # Aggregate per peak:
    # - repeat_class: LINE/SINE/LTR if only one class; otherwise "complex"
    # - repeat_name / repeat_family: comma-separated unique values
    agg <- df_hits[
      ,
      .(
        repeat_class  = if (uniqueN(repClass) == 1) repClass[1] else "complex",
        repeat_name   = paste(unique(repName),   collapse = ","),
        repeat_family = paste(unique(repFamily), collapse = ",")
      ),
      by = q
    ]
    
    # Fill the results back into res (q indexes correspond to row order of res/peak_gr)
    res$repeat_class[agg$q]  <- agg$repeat_class
    res$repeat_name[agg$q]   <- agg$repeat_name
    res$repeat_family[agg$q] <- agg$repeat_family
  }
  # regenerate rownames based on peak coordinates
  res$peak_id <- paste0(res$chr, ":", res$start, "-", res$end)
  rownames(res) <- res$peak_id
  
  
  annotated_results[[mark]] <- res
  cat("✅ Completed:", mark, "\n")
}

saveRDS(annotated_results, file = annotated_rds)

