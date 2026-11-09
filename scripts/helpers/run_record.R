# Run records for diagnostic-dependent parameters (e.g. WGCNA soft power).
# A run record (config/runs/<run_id>.json) holds a value chosen by a person from
# a diagnostic plot, who chose it and when, and a fingerprint of the input the
# choice was made on. Downstream modules refuse to run if the value is missing
# or the fingerprint no longer matches the current input.
# Requires: jsonlite, digest, here

run_record_path <- function(cfg) here::here("config", "runs", paste0(cfg$run_id, ".json"))

is_missing <- function(x) is.null(x) || length(x) == 0 || all(is.na(x))

# Fingerprint of a WGCNA input matrix (samples x genes)
wgcna_fingerprint <- function(datExpr, gene_types) {
  list(
    n_samples  = nrow(datExpr),
    n_genes    = ncol(datExpr),
    gene_types = if (is.null(gene_types)) "all" else paste(sort(gene_types), collapse = ","),
    digest     = digest::digest(list(rownames(datExpr), colnames(datExpr)), algo = "md5")
  )
}

fingerprints_match <- function(a, b) {
  all(vapply(c("n_samples", "n_genes", "gene_types", "digest"),
             function(k) identical(as.character(a[[k]]), as.character(b[[k]])), logical(1)))
}

read_run_record <- function(cfg) {
  path <- run_record_path(cfg)
  if (!file.exists(path)) {
    stop("No run record at ", path, ". Run Module 3A first.", call. = FALSE)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

# Creates a run record with an empty selection. Never overwrites an existing file.
create_run_record <- function(cfg, fingerprint) {
  path <- run_record_path(cfg)
  stopifnot(!file.exists(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  rec <- list(
    run_id = cfg$run_id,
    wgcna_power = list(
      selected_value      = NA,
      selected_by         = NA,
      selected_on         = NA,
      selection_criterion = cfg$wgcna$soft_power$selection_criterion,
      r2_threshold        = cfg$wgcna$soft_power$r2_threshold,
      diagnostic_artifact = "results/WGCNA/SoftThreshold.pdf",
      input_fingerprint   = fingerprint
    )
  )
  jsonlite::write_json(rec, path, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
  path
}

# Stops unless a power has been recorded with attribution
require_selected_power <- function(rec) {
  p <- rec$wgcna_power
  if (is_missing(p$selected_value) || is_missing(p$selected_by) || is_missing(p$selected_on)) {
    stop("Soft-thresholding power not recorded. Inspect results/WGCNA/SoftThreshold.pdf ",
         "from Module 3A, then fill selected_value, selected_by and selected_on in ",
         "config/runs/<run_id>.json.", call. = FALSE)
  }
  if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", p$selected_on) || is.na(as.Date(p$selected_on, optional = TRUE))) {
    stop("selected_on must be a date in YYYY-MM-DD format (got '", p$selected_on, "').", call. = FALSE)
  }
  as.numeric(p$selected_value)
}

# Stops if the input has changed since the power was selected
assert_fingerprint <- function(rec, current) {
  recorded <- rec$wgcna_power$input_fingerprint
  if (!fingerprints_match(recorded, current)) {
    stop("WGCNA input differs from the input the recorded power was selected on.\n",
         "  recorded: ", paste(names(recorded), unlist(recorded), sep = "=", collapse = "; "), "\n",
         "  current:  ", paste(names(current), unlist(current), sep = "=", collapse = "; "), "\n",
         "Re-run Module 3A, re-read the diagnostic plot, and update the run record.",
         call. = FALSE)
  }
  invisible(TRUE)
}
