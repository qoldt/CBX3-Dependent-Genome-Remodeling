# ================================================================
# MINUTE dependency setup - MOVED to ../shared/setup.R
# ----------------------------------------------------------------
# The manifest and installer are shared across assays now (ChIP, GAM, ...), so
# they live in shared/setup.R. This shim keeps the documented entry points
# working unchanged from the ChIP/ directory:
#     Rscript setup.R        install/repair dependencies
#     source("setup.R")      definitions only (what config.R does)
# ================================================================
source("../shared/setup.R")

# shared/setup.R only self-installs when IT is the script being run; sourced
# from here it just defines. So re-check at this level: if THIS file was run
# directly, do the install.
if (sys.nframe() == 0L) {
  minute_install_deps()
  message("All dependencies present:\n  ", paste(minute_all_deps(), collapse = ", "))
}
