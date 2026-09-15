test_that("run_ldx_pipeline applies an explicit multiallelic drop policy", {
  source_vcf <- system.file(
    "extdata", "example_genotypes.vcf", package = "HapBlockR"
  )
  skip_if_not(file.exists(source_vcf), "example VCF is not installed")

  work_dir <- tempfile("hapblockr_multiallelic_pipeline_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  input_vcf <- file.path(work_dir, "multiallelic.vcf")
  lines <- readLines(source_vcf, warn = FALSE)
  data_rows <- which(!startsWith(lines, "#"))
  fields <- strsplit(lines[data_rows[1L]], "\t", fixed = TRUE)[[1L]]
  dropped_id <- fields[3L]
  fields[5L] <- paste(fields[5L], "G", sep = ",")
  lines[data_rows[1L]] <- paste(fields, collapse = "\t")
  writeLines(lines, input_vcf)

  expect_error(
    run_ldx_pipeline(
      geno_source = input_vcf,
      out_dir = work_dir,
      out_blocks = file.path(work_dir, "error_blocks.csv"),
      out_diversity = file.path(work_dir, "error_diversity.csv"),
      out_hap_matrix = file.path(work_dir, "error_haplotypes.tsv"),
      phase = FALSE,
      multiallelic = "error",
      verbose = FALSE
    ),
    "multiallelic"
  )

  result <- run_ldx_pipeline(
    geno_source = input_vcf,
    out_dir = work_dir,
    out_blocks = file.path(work_dir, "blocks.csv"),
    out_diversity = file.path(work_dir, "diversity.csv"),
    out_hap_matrix = file.path(work_dir, "haplotypes.tsv"),
    phase = FALSE,
    multiallelic = "drop",
    maf_cut = 0,
    min_callrate = 0,
    CLQcut = 0.5,
    min_snps_chr = 10L,
    min_snps_block = 3L,
    use_bigmemory = FALSE,
    verbose = FALSE
  )

  expect_false(dropped_id %in% result$snp_info_filtered$SNP)
  expect_equal(nrow(result$snp_info_filtered), length(data_rows) - 1L)
  expect_gt(nrow(result$blocks), 0L)
})
