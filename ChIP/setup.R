# ================================================================
# MINUTE pipeline - dependency manifest + installer
# ----------------------------------------------------------------
# SINGLE SOURCE OF TRUTH for which R packages the pipeline needs.
# config.R sources this file and calls minute_ensure_deps(), so a fresh
# clone installs what it is missing on the first `Rscript run_MINUTE.R`.
#
#   Rscript setup.R          install/repair dependencies, then exit
#   MINUTE_AUTO_INSTALL=0    do NOT install from config.R; missing packages
#                            raise an error listing them instead
#
# Sourcing this file has no side effects beyond defining these objects.
# Adding a library() to config.R? Add the package to the manifest below too.
# ================================================================

# --- Manifest: package -> where it comes from --------------------
# (base/recommended packages such as `parallel` are deliberately absent)
minute_deps <- list(
  # ragg: rasterisation backend for the big ComplexHeatmap bodies - see
  # ht_raster_device in config.R for why the default (cairo) is not usable.
  cran = c("data.table", "dplyr", "ggplot2", "ggrepel", "circlize", "ltc", "ragg"),
  bioc = c("GenomicRanges", "GenomeInfoDb", "rtracklayer", "DESeq2",
           "AnnotationDbi", "ComplexHeatmap", "ChIPseeker",
           "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm39.knownGene"),
  # GitHub packages: "package name" = "user/repo@ref"
  # PINNED - see minute_pins below. Never relax this to a bare "user/repo":
  # that tracks main, and a wigglescout change silently alters every number the
  # pipeline produces (bw_loci() is how all signal is read).
  github = c(wigglescout = "cnluzon/wigglescout@01b0988e010d97107556202f3ef7e952a6826037")
)

# --- Version pins ------------------------------------------------
# Exact versions the pipeline is validated against. A pinned package that is
# installed at a DIFFERENT version counts as missing, so setup.R reinstalls the
# pin rather than silently accepting whatever is on the machine.
#
# wigglescout 0.21.2 == cnluzon/wigglescout main @ 01b0988 (2026-07). The repo
# publishes no git tags, so the ref above is a commit SHA - the only stable
# handle it offers. To move the pin: install the new version, check
# packageVersion("wigglescout") and packageDescription("wigglescout")$RemoteSha,
# update BOTH lines, and re-run the pipeline to confirm the numbers hold.
#
# History: the repo used to pin furrr 0.2.3 / future 1.23.0 / globals 0.14.0
# instead. Those pinned wigglescout's parallel map indirectly - and stopped
# meaning anything when 0.21 dropped future/furrr, which left wigglescout itself
# unpinned. Pin the package that is actually read, not its dependencies.
minute_pins <- c(wigglescout = "0.21.2")

# CRAN/Bioconductor packages are NOT pinned individually: BiocManager already
# ties the Bioconductor release to the R version (validated on R 4.6.1 /
# Bioc 3.23). If you need the full graph locked, add renv - this manifest is a
# floor, not a lockfile.

minute_all_deps <- function() {
  c(minute_deps$cran, minute_deps$bioc, names(minute_deps$github))
}

# Is this dependency satisfied? Installed - and, if pinned, at the pinned
# version. A pinned package sitting at the wrong version is NOT satisfied: that
# is the whole point of the pin, so it gets reinstalled rather than accepted.
minute_dep_ok <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  if (!p %in% names(minute_pins)) return(TRUE)
  identical(as.character(utils::packageVersion(p)), minute_pins[[p]])
}

# Manifest packages that are missing or off-pin (in load order).
minute_missing_deps <- function() Filter(function(p) !minute_dep_ok(p), minute_all_deps())

# Human-readable reason per unsatisfied package, for messages.
minute_dep_problem <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(paste0(p, " (not installed)"))
  paste0(p, " (", utils::packageVersion(p), " installed, pinned to ", minute_pins[[p]], ")")
}

# A library path we can actually write to. R_LIBS_USER is the per-user default
# but does NOT exist until something creates it - a fresh macOS/R install has an
# empty, root-owned system library and no user library at all, which is the most
# common "install fails" case here.
minute_lib <- function() {
  writable <- Filter(function(p) dir.exists(p) && file.access(p, 2) == 0, .libPaths())
  if (length(writable)) return(writable[[1]])
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib) || user_lib == "NULL") {
    user_lib <- file.path(path.expand("~"), "R", paste0(R.version$platform, "-library"),
                          paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = "."))
  }
  user_lib <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]][1]
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(user_lib)) stop("No writable R library found; tried: ", user_lib)
  .libPaths(c(user_lib, .libPaths()))
  message("Created R library: ", user_lib)
  user_lib
}

# Install every missing dependency. Called by minute_ensure_deps(); also usable
# directly to repair an environment (`Rscript setup.R`).
minute_install_deps <- function(lib = minute_lib()) {
  if (!nzchar(getOption("repos")["CRAN"]) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  options(Ncpus = max(1L, parallel::detectCores() - 1L))
  need <- function(p) !requireNamespace(p, quietly = TRUE)

  for (p in c("BiocManager", "remotes")) if (need(p)) install.packages(p, lib = lib)

  todo <- Filter(need, minute_deps$cran)
  if (length(todo)) {
    message("Installing from CRAN: ", paste(todo, collapse = ", "))
    install.packages(todo, lib = lib)
  }

  todo <- Filter(need, minute_deps$bioc)
  if (length(todo)) {
    message("Installing from Bioconductor ", BiocManager::version(), ": ",
            paste(todo, collapse = ", "))
    BiocManager::install(todo, lib = lib, ask = FALSE, update = FALSE)
  }

  # force = TRUE: an off-pin copy is already installed, so remotes would other-
  # wise consider the requirement met and skip it.
  todo <- names(minute_deps$github)[!vapply(names(minute_deps$github), minute_dep_ok, logical(1))]
  for (p in todo) {
    message("Installing from GitHub: ", minute_deps$github[[p]])
    remotes::install_github(minute_deps$github[[p]], lib = lib,
                            upgrade = "never", force = TRUE)
  }

  still <- minute_missing_deps()
  if (length(still)) {
    stop("Failed to install: ", paste(vapply(still, minute_dep_problem, ""), collapse = ", "),
         "\nInstall them by hand, then re-run. Library used: ", lib)
  }
  invisible(TRUE)
}

# Called from config.R on every run. Fast no-op once everything is present:
# requireNamespace() on installed packages is cheap and does not attach them.
minute_ensure_deps <- function() {
  missing <- minute_missing_deps()
  if (!length(missing)) return(invisible(TRUE))
  detail <- paste(vapply(missing, minute_dep_problem, ""), collapse = ", ")
  if (identical(Sys.getenv("MINUTE_AUTO_INSTALL", "1"), "0")) {
    stop("Unsatisfied R packages: ", detail,
         "\nRun `Rscript setup.R` (or unset MINUTE_AUTO_INSTALL) to install them.")
  }
  message("Unsatisfied R packages: ", detail)
  message("Installing them now (MINUTE_AUTO_INSTALL=0 disables this)...")
  minute_install_deps()
}

# `Rscript setup.R` -> install; `source("setup.R")` -> definitions only.
if (sys.nframe() == 0L) {
  minute_install_deps()
  message("All MINUTE dependencies present:\n  ",
          paste(minute_all_deps(), collapse = ", "))
}
