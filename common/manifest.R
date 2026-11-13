# Run manifests (md5 + dims) and run-context logging for every result directory.

# Paths are stored relative to dirname(out_file) when under it, else absolute.
write_manifest <- function(paths, out_file, role = c("input", "output")) {
  role  <- match.arg(role)
  files <- unlist(lapply(paths, function(p)
    if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE) else p))
  files <- normalizePath(files, mustWork = TRUE)
  files <- files[files != normalizePath(out_file, mustWork = FALSE)]
  root  <- normalizePath(dirname(out_file), mustWork = FALSE)

  dims <- t(vapply(files, function(f) {
    obj <- tryCatch(
      if (grepl("\\.rds$", f, ignore.case = TRUE)) readRDS(f)
      else if (grepl("\\.csv$", f, ignore.case = TRUE)) as.data.frame(data.table::fread(f))
      else NULL,
      error = function(e) NULL)
    if (is.data.frame(obj) || is.matrix(obj)) c(nrow(obj), ncol(obj)) else c(NA, NA)
  }, numeric(2)))

  info <- file.info(files)
  rel  <- ifelse(startsWith(files, paste0(root, "/")), substring(files, nchar(root) + 2), files)
  m <- data.frame(role = role, path = rel, bytes = info$size,
                  md5 = unname(tools::md5sum(files)),
                  mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
                  n_rows = dims[, 1], n_cols = dims[, 2],
                  stringsAsFactors = FALSE, row.names = NULL)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  write.table(m, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(m)
}

verify_manifest <- function(manifest_tsv, base_dir) {
  m    <- read.delim(manifest_tsv, stringsAsFactors = FALSE)
  full <- ifelse(startsWith(m$path, "/"), m$path, file.path(base_dir, m$path))
  ok   <- file.exists(full)
  ok[ok] <- unname(tools::md5sum(full[ok])) == m$md5[ok]
  if (all(ok)) cat(sprintf("OK %d/%d\n", sum(ok), length(ok)))
  else cat(sprintf("MISMATCH %d/%d:\n", sum(!ok), length(ok)), paste0("  ", m$path[!ok], "\n"), sep = "")
  invisible(all(ok))
}

log_run_context <- function(config, out_dir) {
  rev_file <- file.path(Sys.getenv("CODE", Sys.getenv("MHSEL_ROOT", ".")), "REVISION")
  pkgs <- c("masbayes", "masreml", "AlphaSimR")
  vers <- vapply(pkgs, function(p)
    tryCatch(as.character(packageVersion(p)), error = function(e) "not installed"), "")
  lines <- c(
    sprintf("time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("host: %s", Sys.info()[["nodename"]]),
    sprintf("SLURM_JOB_ID: %s | SLURM_ARRAY_TASK_ID: %s",
            Sys.getenv("SLURM_JOB_ID", "-"), Sys.getenv("SLURM_ARRAY_TASK_ID", "-")),
    sprintf("REVISION: %s", if (file.exists(rev_file)) readLines(rev_file, 1) else "unknown"),
    sprintf("packages: %s", paste(pkgs, vers, collapse = " | ")),
    "config:", deparse(config),
    "sessionInfo:", capture.output(sessionInfo()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "run_context.txt"))
  invisible(lines)
}
