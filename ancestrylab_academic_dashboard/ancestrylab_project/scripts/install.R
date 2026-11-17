options(repos = c(CRAN = "https://cloud.r-project.org"))
required <- c("remotes", "shiny", "bslib", "DT", "yaml", "jsonlite", "ggplot2", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

# Pin ADMIXTOOLS 2 to a reviewed upstream commit for reproducible installs.
# Set ADMIXTOOLS_REF to another reviewed commit only when intentionally upgrading.
admixtools_ref <- Sys.getenv("ADMIXTOOLS_REF", unset = "f3969126a0ceca3f05d1c9a11c10f10df2b7024e")
if (!requireNamespace("admixtools", quietly = TRUE)) {
  remotes::install_github(paste0("uqrmaie1/admixtools@", admixtools_ref), upgrade = "never")
}
message("AncestryLab dependencies installed. ADMIXTOOLS 2: ", as.character(utils::packageVersion("admixtools")))
source("R/engine.R")
stopifnot(isTRUE(tryCatch({ validate_population_design("T", c("A", "B"), c("R1", "R2")); TRUE }, error=function(e) FALSE)))
