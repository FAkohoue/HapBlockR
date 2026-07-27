#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE))
  stop("The benchmark requires devtools.")
if (!requireNamespace("digest", quietly = TRUE))
  stop("The benchmark requires digest.")

devtools::load_all(quiet = TRUE)
set.seed(20260725)
n_samples <- 300L
n_variants <- 600L
n_repetitions <- 3L

geno <- matrix(
  stats::rbinom(n_samples * n_variants, 2L, 0.35),
  nrow = n_samples,
  ncol = n_variants,
  dimnames = list(
    sprintf("G%04d", seq_len(n_samples)),
    sprintf("S%05d", seq_len(n_variants))
  )
)
missing <- sample.int(length(geno), floor(0.01 * length(geno)))
geno[missing] <- NA_real_

elapsed <- function(threads) {
  replicate(
    n_repetitions,
    unname(system.time(
      result <- HapBlockR::compute_r2(
        geno,
        n_threads = threads
      )
    )[["elapsed"]])
  )
}

time_one <- elapsed(1L)
time_two <- elapsed(2L)
result_one <- HapBlockR::compute_r2(geno, n_threads = 1L)
result_two <- HapBlockR::compute_r2(geno, n_threads = 2L)

reference_geno <- geno[, seq_len(30L), drop = FALSE]
for (column in seq_len(ncol(reference_geno))) {
  missing_column <- is.na(reference_geno[, column])
  if (any(missing_column))
    reference_geno[missing_column, column] <- mean(
      reference_geno[, column], na.rm = TRUE
    )
}
reference <- stats::cor(reference_geno)^2
diag(reference) <- 0
compiled_reference <- result_one[seq_len(30L), seq_len(30L)]

peak_rss_mb <- NA_real_
status_path <- "/proc/self/status"
if (file.exists(status_path)) {
  status <- readLines(status_path, warn = FALSE)
  vmhwm <- grep("^VmHWM:", status, value = TRUE)
  if (length(vmhwm)) {
    peak_rss_mb <- as.numeric(
      sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", vmhwm[1L])
    ) / 1024
  }
}

compiler_command <- if (.Platform$OS.type == "windows") {
  file.path(R.home("bin"), "x64", "Rcmd.exe")
} else {
  file.path(R.home("bin"), "R")
}
compiler_arguments <- if (.Platform$OS.type == "windows") {
  c("config", "CXX17")
} else {
  c("CMD", "config", "CXX17")
}
compiler <- tryCatch(
  suppressWarnings(system2(
    compiler_command,
    compiler_arguments,
    stdout = TRUE,
    stderr = TRUE
  )),
  error = function(e) NA_character_
)
if (
  !length(compiler) ||
  all(is.na(compiler)) ||
  !is.null(attr(compiler, "status")) && attr(compiler, "status") != 0L
) {
  compiler_executable <- Sys.which(c("g++", "clang++"))
  compiler_executable <- compiler_executable[nzchar(compiler_executable)][1L]
  compiler <- if (length(compiler_executable)) {
    tryCatch(
      suppressWarnings(system2(
        compiler_executable,
        "--version",
        stdout = TRUE,
        stderr = TRUE
      ))[1L],
      error = function(e) NA_character_
    )
  } else {
    NA_character_
  }
}
processor <- if (file.exists("/proc/cpuinfo")) {
  cpu <- grep(
    "^model name[[:space:]]*:",
    readLines("/proc/cpuinfo", warn = FALSE),
    value = TRUE
  )
  if (length(cpu)) sub("^.*:[[:space:]]*", "", cpu[1L]) else NA_character_
} else {
  Sys.info()[["machine"]]
}

output <- data.frame(
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_sha256 = digest::digest(geno, algo = "sha256"),
  n_samples = n_samples,
  n_variants = n_variants,
  missing_rate = mean(is.na(geno)),
  repetitions = n_repetitions,
  median_elapsed_one_thread_seconds = stats::median(time_one),
  median_elapsed_two_threads_seconds = stats::median(time_two),
  elapsed_one_thread_sd = stats::sd(time_one),
  elapsed_two_threads_sd = stats::sd(time_two),
  peak_resident_mb = peak_rss_mb,
  one_two_thread_max_abs_error = max(abs(result_one - result_two)),
  scalar_reference_max_abs_error =
    max(abs(compiled_reference - reference)),
  correctness_passed =
    isTRUE(all.equal(result_one, result_two, tolerance = 1e-12)) &&
    isTRUE(all.equal(
      unname(compiled_reference),
      unname(reference),
      tolerance = 1e-10
    )),
  R_version = R.version.string,
  HapBlockR_version = as.character(utils::packageVersion("HapBlockR")),
  compiler = paste(compiler, collapse = " "),
  operating_system = paste(Sys.info()[c("sysname", "release")],
                           collapse = " "),
  processor = processor,
  stringsAsFactors = FALSE
)

thresholds <- utils::read.csv(
  "benchmarks/thresholds.csv",
  stringsAsFactors = FALSE
)
for (index in seq_len(nrow(thresholds))) {
  observed <- output[[thresholds$metric[index]]]
  if (is.finite(observed) && observed > thresholds$maximum[index])
    stop(
      thresholds$metric[index], " = ", observed,
      " exceeds ", thresholds$maximum[index], "."
    )
}
if (!output$correctness_passed)
  stop("Compiled LD output failed the correctness comparison.")

dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  output,
  "benchmarks/results/current.csv",
  row.names = FALSE,
  na = ""
)
print(output)
