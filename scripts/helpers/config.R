# Shared configuration, path and logging helpers. Sourced by every module.
# Requires: jsonlite, here, knitr

.log_env <- new.env(parent = emptyenv())

# ---- Configuration ----------------------------------------------------------

# Reads config/config.json, then overlays config/config_local.json if present.
# Creates the standard project_dir layout.
load_config <- function(config_dir = here::here("config")) {
  cfg <- jsonlite::fromJSON(file.path(config_dir, "config.json"))
  local_file <- file.path(config_dir, "config_local.json")
  if (file.exists(local_file)) {
    cfg <- utils::modifyList(cfg, jsonlite::fromJSON(local_file))
  }

  project_dir <- cfg$paths$project_dir
  if (is.null(project_dir) || grepl("^/path/to", project_dir)) {
    stop("paths$project_dir is not set. Create config/config_local.json containing:\n",
         '  {"paths": {"project_dir": "/absolute/path/to/project"}}', call. = FALSE)
  }
  cfg$paths$project_dir <- normalizePath(path.expand(project_dir), mustWork = FALSE)

  for (d in c("data/raw", "data/external", "data/processed", "results", "logs")) {
    dir.create(file.path(cfg$paths$project_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
  cfg
}

# Builds a path inside project_dir.
proj_path <- function(cfg, ...) file.path(cfg$paths$project_dir, ...)

# Reports files that are being reused instead of regenerated.
note_existing <- function(paths) {
  for (p in paths) {
    message("Existing file used, not regenerated: ", p,
            " (modified ", format(file.mtime(p), "%Y-%m-%d %H:%M"), ")")
  }
}

# ---- Logging ----------------------------------------------------------------

git_state <- function() {
  run_git <- function(args) {
    tryCatch(suppressWarnings(system2("git", c("-C", shQuote(here::here()), args),
                                      stdout = TRUE, stderr = FALSE)),
             error = function(e) character(0))
  }
  commit <- run_git(c("rev-parse", "--short", "HEAD"))
  list(commit = if (length(commit) == 1) commit else "none (not a git commit)",
       dirty  = length(run_git(c("status", "--porcelain"))) > 0)
}

# Opens a timestamped log in <project_dir>/logs/.
#   Interactive runs: console output is split to the log with sink();
#                     messages and warnings are appended by global handlers.
#   Knitting:         knitr hooks append source, output, messages, warnings
#                     and errors from all subsequent chunks to the log.
start_log <- function(cfg, module) {
  log_path <- proj_path(cfg, "logs",
                        sprintf("%s_%s.log", module, format(Sys.time(), "%Y%m%d_%H%M%S")))
  knitting <- isTRUE(getOption("knitr.in.progress"))
  git <- git_state()

  writeLines(c(
    sprintf("Module:      %s", module),
    sprintf("Started:     %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Project dir: %s", cfg$paths$project_dir),
    sprintf("Git commit:  %s%s", git$commit,
            if (git$dirty) " (uncommitted changes present)" else ""),
    sprintf("Mode:        %s", if (knitting) "knit" else "interactive"),
    strrep("=", 72)
  ), log_path)

  .log_env$path <- log_path
  .log_env$knitting <- knitting
  .log_env$active <- TRUE

  if (knitting) {
    for (hook in c("source", "output", "message", "warning", "error")) {
      local({
        h <- hook
        default_hook <- knitr::knit_hooks$get(h)
        knitr::knit_hooks$set(stats::setNames(list(function(x, options) {
          if (isTRUE(.log_env$active)) {
            text <- if (h == "source") paste0("> ", x) else x
            cat(text, sep = "\n", file = .log_env$path, append = TRUE)
          }
          default_hook(x, options)
        }), h))
      })
    }
  } else {
    .log_env$con <- file(log_path, open = "a")
    sink(.log_env$con, split = TRUE)
    tryCatch(
      globalCallingHandlers(
        message = function(m) cat("[message] ", conditionMessage(m), file = .log_env$con),
        warning = function(w) cat("[warning] ", conditionMessage(w), "\n", file = .log_env$con)
      ),
      error = function(e) message("Messages/warnings not captured in log: ", conditionMessage(e))
    )
  }
  message("Logging to ", log_path)
  invisible(log_path)
}

# Writes sessionInfo() to the log and closes it.
stop_log <- function() {
  if (!isTRUE(.log_env$active)) return(invisible(NULL))
  if (.log_env$knitting) {
    cat(c(strrep("=", 72), utils::capture.output(utils::sessionInfo())),
        sep = "\n", file = .log_env$path, append = TRUE)
    .log_env$active <- FALSE
    print(utils::sessionInfo())
  } else {
    cat(strrep("=", 72), "\n")
    print(utils::sessionInfo())
    globalCallingHandlers(NULL)
    while (sink.number() > 0) sink()
    close(.log_env$con)
    .log_env$active <- FALSE
  }
  message("Log closed: ", .log_env$path)
  invisible(.log_env$path)
}
