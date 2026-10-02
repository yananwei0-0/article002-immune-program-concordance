sci27_repo_root <- normalizePath(Sys.getenv("SCI27_REPO_ROOT", getwd()), mustWork = FALSE)

sci27_legacy <- function(...) file.path(.Platform$file.sep, ...)

sci27_expand_legacy_paths <- function(value) {
  external <- file.path(sci27_repo_root, "external_data")
  mapping <- c(
    SCI27_RESEARCH_ROOT = Sys.getenv("SCI27_RESEARCH_ROOT", sci27_repo_root),
    SCI27_DATA_ROOT = Sys.getenv("SCI27_DATA_ROOT", external),
    SCI27_BIOBANK_ROOT = Sys.getenv("SCI27_BIOBANK_ROOT", file.path(external, "biobank")),
    SCI27_TMP_ROOT = Sys.getenv("SCI27_TMP_ROOT", tempdir())
  )
  legacy <- c(
    SCI27_RESEARCH_ROOT = sci27_legacy("Users", "yanan", "Desktop", "research"),
    SCI27_DATA_ROOT = sci27_legacy("Users", "yanan", "Desktop", "bioinformatics_analysis"),
    SCI27_BIOBANK_ROOT = sci27_legacy("Volumes", "Biobank"),
    SCI27_TMP_ROOT = sci27_legacy("private", "tmp")
  )
  legacy <- c(legacy, SCI27_TMP_ROOT_SHORT = sci27_legacy("tmp"))
  replacements <- c(mapping, SCI27_TMP_ROOT_SHORT = mapping[["SCI27_TMP_ROOT"]])
  tokens <- paste0("__SCI27_PATH_", seq_along(legacy), "__")
  result <- value
  for (index in seq_along(legacy)) result <- gsub(legacy[[index]], tokens[[index]], result, fixed = TRUE)
  for (index in seq_along(tokens)) result <- gsub(tokens[[index]], replacements[[index]], result, fixed = TRUE)
  result
}
