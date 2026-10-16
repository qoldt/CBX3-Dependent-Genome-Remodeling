# ================================================================
# GAM_7 - Local contact counts from the PAIRS, and distance-stratified coupling
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_7_Pairs_Local_Contacts.R     (after GAM_1)
#
# WHY THIS EXISTS. Every earlier stage used the delivered Hotspots track, whose
# per-bin value is the number of significant differential contacts within
# 0-5 Mb OF THAT BIN (it is in the filename: dist0-5000000). That value is not a
# local measurement - it integrates a 5 Mb neighbourhood. GAM_3 therefore
# correlated a genuinely 20 kb variable (ChIP log2FC) against one smeared over
# 5 Mb, a ~250x scale mismatch, and reported a null. The null is real for the
# question as posed, but it cannot distinguish "no relationship" from "the
# contact variable cannot resolve the scale the chromatin variable varies on".
#
# This stage rebuilds the contact signal from Diff_contacts/, where each pair
# carries both endpoints and its separation. A bin's count here is the number of
# differential contacts with an ENDPOINT IN THAT BIN - genuinely local - and can
# be split by distance, which the Hotspots summary cannot.
#
# Memory: ~25M pairs per direction genome-wide, so everything is per chromosome
# and only the needed columns are read.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
  "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(g1)
a <- as.data.table(x$annot); setorder(a, chr, start)
for (mk in colnames(x$l2fc)) a[[paste0("l2fc_", mk)]] <- x$l2fc[match(a$bin_id, rownames(x$l2fc)), mk]

# Distance bands. The Hotspots track collapses everything under 5 Mb into one
# number and discards everything above it; the pairs run to 57 Mb. These bands
# separate the scales a compaction change and a compartment change would act on.
BANDS <- list("0.02-1Mb" = c(2e4, 1e6), "1-5Mb" = c(1e6, 5e6),
              "5-10Mb" = c(5e6, 1e7), ">10Mb" = c(1e7, Inf))

bin_index <- function(chrom, pos) {
  # pos is a 0-based start from the bg2; bins are gam_binsize wide from 0
  paste0(chrom, ":", (pos %/% gam_binsize) * gam_binsize)
}
a[, key := paste0(chr, ":", start)]
setkey(a, key)

count_one <- function(chrom, dir) {
  f <- diff_file(chrom, dir)
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = c("chrom_x","start_x","chrom_y","start_y","bindist"),
             showProgress = FALSE)
  if (!nrow(d)) return(NULL)
  d[, band := cut(bindist, breaks = c(0, 1e6, 5e6, 1e7, Inf),
                  labels = names(BANDS), right = FALSE)]
  # A pair contributes to BOTH of its endpoint bins: the question is "how many
  # differential contacts does this bin participate in", so each end counts.
  both <- rbindlist(list(
    d[, .(key = bin_index(chrom_x, start_x), band)],
    d[, .(key = bin_index(chrom_y, start_y), band)]))
  out <- both[, .N, by = .(key, band)]
  rm(d, both); gc(FALSE)
  out
}

cat("Counting pairs per bin, per distance band (", length(gam_autosomes), "chromosomes )\n")
res <- list()
for (ch in gam_autosomes) {
  for (dir in c("lost", "gained")) {
    o <- count_one(ch, dir)
    if (!is.null(o)) { o[, dirn := dir]; res[[paste(ch, dir)]] <- o }
  }
  cat("  ", ch, "done\n")
}
cnt <- rbindlist(res); rm(res); gc(FALSE)

# wide: one column per direction x band
w <- dcast(cnt, key ~ dirn + band, value.var = "N", fill = 0)
a <- merge(a, w, by = "key", all.x = TRUE)
for (j in setdiff(names(w), "key")) set(a, which(is.na(a[[j]])), j, 0L)

a[, local_lost   := rowSums(.SD), .SDcols = grep("^lost_",   names(a), value = TRUE)]
a[, local_gained := rowSums(.SD), .SDcols = grep("^gained_", names(a), value = TRUE)]

cat("\n=== local (pairs) vs delivered (0-5Mb window) contact counts ===\n")
cat(sprintf("  delivered lost : median %4d  IQR %5.0f\n", median(a$n_contacts_lost), IQR(a$n_contacts_lost)))
cat(sprintf("  local     lost : median %4d  IQR %5.0f\n", median(a$local_lost), IQR(a$local_lost)))
cat(sprintf("  rho(delivered, local) = %+.3f  <- how much the 5 Mb window smears\n",
            cor(a$n_contacts_lost, a$local_lost, method = "spearman")))

# --- the re-test: does a LOCAL contact measure couple to the ChIP change? -----
marks <- colnames(x$l2fc)
cp <- rbindlist(lapply(marks, function(mk) {
  v <- a[[paste0("l2fc_", mk)]]
  rbindlist(lapply(c(list(all = NULL), BANDS), function(bd) NULL))  # placeholder
}))
rows <- list()
for (mk in marks) {
  v <- a[[paste0("l2fc_", mk)]]
  for (dir in c("lost","gained")) {
    rows[[length(rows)+1]] <- data.table(mark = mk, direction = dir, band = "all (local)",
      rho = cor(v, a[[paste0("local_", dir)]], method = "spearman", use = "complete.obs"))
    for (bn in names(BANDS)) {
      col <- paste0(dir, "_", bn)
      if (col %in% names(a))
        rows[[length(rows)+1]] <- data.table(mark = mk, direction = dir, band = bn,
          rho = cor(v, a[[col]], method = "spearman", use = "complete.obs"))
    }
  }
  rows[[length(rows)+1]] <- data.table(mark = mk, direction = "delivered", band = "0-5Mb window",
    rho = cor(v, a$n_contacts_lost, method = "spearman", use = "complete.obs"))
}
cp <- rbindlist(rows)
fwrite(cp, file.path(tables_dir, "gam7_coupling_by_distance.tsv"), sep = "\t")

cat("\n=== coupling of per-bin ChIP log2FC to LOCAL contact counts, by distance ===\n")
pr <- dcast(cp[direction != "delivered"], mark + direction ~ band, value.var = "rho")
print(pr, row.names = FALSE, digits = 3)
cat("\n  (GAM_3's delivered-track value, for comparison)\n")
print(cp[direction == "delivered", .(mark, rho = round(rho, 3))], row.names = FALSE)

saveRDS(list(bins = a, coupling = cp, bands = BANDS),
        file.path(rds_dir, "gam7_local_contacts.rds"))
fwrite(a[, c("bin_id","chr","start","end","local_lost","local_gained",
             grep("^(lost|gained)_", names(a), value = TRUE)), with = FALSE],
       file.path(tables_dir, "gam7_local_contact_counts.tsv.gz"), sep = "\t")
cat("\nGAM_7 complete.\n")
