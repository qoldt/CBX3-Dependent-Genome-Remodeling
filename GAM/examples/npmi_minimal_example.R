# ================================================================
# npmi_minimal_example.R - what GAM_2 does to an NPMI matrix, one step at a time
# ----------------------------------------------------------------
#   cd GAM && Rscript examples/npmi_minimal_example.R     (~20 s, ~1 GB RAM)
#
# Teaching script, not a pipeline stage. ONE genotype (WT), ONE chromosome
# (chr19, the smallest), ONE bin - no loops. GAM_2_NPMI_Features.R does the same
# thing for every bin x chromosome x genotype. The last section checks the
# numbers below against GAM_2's own cached output for this bin.
# Only needs data.table; config.R / utils.R are not sourced so every piece is visible.
# ================================================================
suppressPackageStartupMessages(library(data.table))

# ---- 1. Where the file is -------------------------------------------------
# utils.R:npmi_file() builds this. Genotype is ONLY in the folder name
# (WT_at20Kb / cKO_at20Kb); the file names are identical between genotypes.
binsize <- 20000L
f <- file.path("NPMI_matrices", "WT_at20Kb",
  "HP1g_mm39.segregation_at20000.2026-02-23.passed_qc.curated.sigma.2.6.NPMI.chr19.squared.txt.gz")
stopifnot(file.exists(f))

# ---- 2. Read it into a square matrix --------------------------------------
# utils.R:load_npmi(). The file is a dense square table: first column = bin
# id, header = the same bin ids. Cell [i, j] = NPMI between bins i and j
# (how much more often they land in the same nuclear slice than chance;
# -1..1, 0 = independent). NA = diagonal or an unmappable bin.
dt  <- fread(f, header = TRUE)
ids <- as.character(dt[[1]])                 # "chr19:0-20000", ...
m   <- as.matrix(dt[, -1, with = FALSE])
rownames(m) <- ids
rm(dt)
n <- nrow(m)
stopifnot(n == ncol(m), identical(colnames(m), ids), isSymmetric(m))
cat(sprintf("chr19 WT: %d x %d bins (%.1f Mb), %.1f%% of cells NA\n",
            n, n, n * binsize / 1e6, 100 * mean(is.na(m))))

# Bins with no data at all (repeats, gaps) are all-NA rows. They stay NA:
# filling them with 0 would mean "never co-segregates", which is a claim.
n_valid <- colSums(!is.na(m))
cat(sprintf("all-NA (unmappable) bins: %d\n", sum(n_valid == 0)))

# ---- 3. Pick one bin ------------------------------------------------------
i <- match("chr19:20000000-20020000", ids)
stopifnot(!is.na(i), n_valid[i] > 0)
v <- m[, i]                                  # this bin's contacts with every other bin
d <- abs(seq_len(n) - i)                     # distance from bin i, in bins

# ---- 4. Short- and long-range contact for that bin ------------------------
# SHORT = partners 1-50 bins away (20 kb - 1 Mb; TAD scale).
# LONG  = partners >= 250 bins away (>= 5 Mb; beyond any TAD).
npmi_short <- mean(v[d >= 1   & d <= 50], na.rm = TRUE)
npmi_long  <- mean(v[d >= 250],           na.rm = TRUE)

# ---- 5. Decay slope --------------------------------------------------------
# How fast contact falls off with distance: mean NPMI in log-spaced distance
# shells (40 kb - 10 Mb), then a straight-line fit of log10(NPMI) on
# log10(distance). Steeper (more negative) = contacts more local.
edges <- unique(round(2 * 1.5^(0:20))); edges <- c(edges[edges < 500], 500)
shell <- cut(d, breaks = edges, right = FALSE)               # [a, b) in bins
shell_mean <- tapply(v, shell, mean, na.rm = TRUE)
shell_n    <- tapply(!is.na(v), shell, sum)
shell_dist <- (head(edges, -1) + tail(edges, -1) - 1) / 2 * binsize  # shell midpoint, bp
use <- !is.na(shell_mean) & shell_mean > 0 & shell_n >= 3     # log needs > 0
decay_slope <- if (sum(use) >= 5)
  unname(coef(lm(log10(shell_mean[use]) ~ log10(shell_dist[use])))[2]) else NA

# ---- 6. Insulation ---------------------------------------------------------
# Mean NPMI in the square that spans bin i: the w bins upstream x the w bins
# downstream. Low = few contacts cross this bin = it acts as a boundary.
# w = 10 bins = 200 kb, the primary window. GAM_2 then subtracts the
# chromosome median (absolute NPMI level differs between datasets; the
# local dip does not), which needs every bin, so the raw value is shown here.
w   <- 10L
blk <- m[(i - w):(i - 1L), (i + 1L):(i + w)]
insulation_raw <- mean(blk, na.rm = TRUE)

# ---- 7. Contact-decay curve for the whole chromosome -----------------------
# Mean NPMI on each diagonal (all bin pairs exactly k bins apart). This is the
# genome-wide version of step 5 and is the curve GAM_2 plots per genotype.
k <- (col(m) - row(m))[upper.tri(m)]
decay_curve <- tapply(m[upper.tri(m)], k, mean, na.rm = TRUE)
cat(sprintf("mean NPMI at 20 kb: %.3f | 1 Mb: %.3f | 10 Mb: %.3f\n",
            decay_curve["1"], decay_curve["50"], decay_curve["500"]))

# ---- 8. The result for this bin -------------------------------------------
res <- data.table(bin_id = ids[i], npmi_short, npmi_long, decay_slope,
                  insulation_raw, n_valid = n_valid[i])
print(res)

# ---- 9. Check against GAM_2's cached reduction ------------------------------
cache <- file.path("results", "rds", "gam2_cache", "chr19_WT.rds")
if (file.exists(cache)) {
  x <- readRDS(cache)
  g <- (if (is.data.frame(x)) x else x$feat)[bin_id == ids[i]]
  stopifnot(isTRUE(all.equal(g$npmi_short,  npmi_short)),
            isTRUE(all.equal(g$npmi_long,   npmi_long)),
            isTRUE(all.equal(g$decay_slope, decay_slope)),
            g$n_valid == n_valid[i])
  cat("matches GAM_2's cached features for this bin\n")
}
