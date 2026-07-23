## tests/testthat/helper.R
## Shared helpers loaded automatically by testthat before every test file.
## -----------------------------------------------------------------------------

#' Skip a test if SimpleMating is installed but too old to export `fn`
#'
#' SimpleMating rewrote its API between 0.1.x (buildCrosses()/selectCrosses()/
#' evoluteCandidates()/etc., no GOCS()/getUsefA()) and 0.2.x (adds GOCS()/
#' getUsefA()/getUsefAD()/getMPV()/getTGV(), which R/ocs.R and
#' R/genomic_mating.R target -- confirmed against the live GitHub main branch
#' NAMESPACE/DESCRIPTION, currently 0.2.1). `skip_if_not_installed("SimpleMating")`
#' alone is not enough to guarantee these tests are safe to run: an installed-
#' but-outdated copy passes that check yet still fails with a real (if now
#' clearly-worded, see R/ocs.R's/R/genomic_mating.R's explicit version-gap
#' check) error. Use alongside `skip_if_not_installed("SimpleMating")`, not
#' instead of it.
#' @param fn Character. Function name expected to be exported (e.g. "GOCS"
#'   for the OCS engine, "getUsefA" for usefulness_criterion()'s
#'   variance_model = "simplemating").
skip_if_simplemating_too_old <- function(fn) {
  if (requireNamespace("SimpleMating", quietly = TRUE) &&
      !(fn %in% getNamespaceExports("SimpleMating"))) {
    testthat::skip(paste0(
      "SimpleMating is installed (version ",
      as.character(utils::packageVersion("SimpleMating")), ") but does not ",
      "export ", fn, "() -- an older, pre-0.2.x release. Reinstall via ",
      "remotes::install_github('Resende-Lab/SimpleMating', force = TRUE) ",
      "to run this test."
    ))
  }
}

#' Skip a test if optiSel is installed but does not export `fn`
#'
#' select_parents_ocs(engine = "optisel") (the true-OCS engine, via
#' R/ocs.R's .run_optisel_ocs()) targets optiSel::candes()/opticont()/
#' matings()/noffspring() directly against that package's own documented
#' API. Use this guard so an installed-but-incompatible optiSel version
#' (missing one of these exports, e.g. a much older release) skips cleanly
#' with an actionable message instead of failing on an unrelated "not an
#' exported object" error.
#' @param fn Character. Function name expected to be exported (one of
#'   "candes", "opticont", "matings", "noffspring").
skip_if_optisel_missing_fn <- function(fn) {
  if (requireNamespace("optiSel", quietly = TRUE) &&
      !(fn %in% getNamespaceExports("optiSel"))) {
    testthat::skip(paste0(
      "optiSel is installed (version ",
      as.character(utils::packageVersion("optiSel")), ") but does not ",
      "export ", fn, "() -- select_parents_ocs(engine = 'optisel') targets ",
      "this function directly. Check ?optiSel::", fn, " against your ",
      "installed version, or update optiSel."
    ))
  }
}

# -- Genotype matrix -----------------------------------------------------------

#' Simulate a random dosage matrix
#' @param n  Number of individuals. Default 60.
#' @param p  Number of SNPs. Default 30.
#' @param seed RNG seed for reproducibility.
#' @return Numeric matrix (n x p), values in {0,1,2}, named rows and columns.
make_geno <- function(n = 60, p = 30, seed = 1L) {
  set.seed(seed)
  G <- matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
  rownames(G) <- paste0("ind", seq_len(n))
  colnames(G) <- paste0("rs",  seq_len(p))
  G
}

# -- SNP metadata --------------------------------------------------------------

#' Build a minimal SNP information data frame
#' @param p    Number of SNPs. Default 30.
#' @param chr  Chromosome label. Default "1".
#' @return Data frame with columns SNP, CHR, POS, REF, ALT. REF and ALT
#'   alternate between A/T and G/C pairs so that nucleotide decoding produces
#'   meaningful output in write_haplotype_character() tests.
make_snpinfo <- function(p = 30, chr = "1") {
  refs <- rep(c("A", "G", "C", "T"), length.out = p)
  alts <- rep(c("T", "C", "A", "G"), length.out = p)
  data.frame(
    SNP = paste0("rs",  seq_len(p)),
    CHR = chr,
    POS = seq(1000L, by = 5000L, length.out = p),
    REF = refs,
    ALT = alts,
    stringsAsFactors = FALSE
  )
}

# -- Block table ---------------------------------------------------------------

#' Build a minimal block table from snp_info
#' @param snp_info  Data frame from make_snpinfo() or ldx_snp_info.
#' @param n_blocks  Number of equal-sized blocks to create. Default 3.
#' @return Data frame compatible with run_Big_LD_all_chr() output: columns
#'   start, end, start.rsID, end.rsID, start.bp, end.bp, CHR, length_bp.
#'   Blocks are equally sized; any remaining SNPs go into the last block.
make_blocks <- function(snp_info, n_blocks = 3L) {
  p        <- nrow(snp_info)
  bsize    <- floor(p / n_blocks)
  starts   <- seq(1L, by = bsize, length.out = n_blocks)
  ends     <- c(starts[-1L] - 1L, p)
  chr      <- as.character(snp_info$CHR[1L])
  rows <- lapply(seq_len(n_blocks), function(i) {
    s <- starts[i]; e <- ends[i]
    data.frame(
      start      = s,
      end        = e,
      start.rsID = snp_info$SNP[s],
      end.rsID   = snp_info$SNP[e],
      start.bp   = snp_info$POS[s],
      end.bp     = snp_info$POS[e],
      CHR        = chr,
      length_bp  = snp_info$POS[e] - snp_info$POS[s] + 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# -- Phenotype / BLUEs ---------------------------------------------------------

#' Simulate pre-adjusted phenotype values (BLUEs) for run_haplotype_prediction()
#' @param geno_matrix  Genotype matrix (individuals x SNPs) from make_geno().
#' @param n_qtl        Number of causal SNPs. Default 5.
#' @param h2           Heritability for signal-to-noise scaling. Default 0.5.
#' @param format       Return format: "vector" (named numeric, default) or
#'   "data.frame" (data frame with columns "id" and "blue").
#' @param seed         RNG seed. Default 42L.
#' @return Named numeric vector or data frame of simulated BLUEs. These
#'   represent the output of a field trial mixed model - the typical input to
#'   GWAS and run_haplotype_prediction().
make_blues <- function(geno_matrix, n_qtl = 5L, h2 = 0.5,
                       format = c("vector", "data.frame"), seed = 42L) {
  format <- match.arg(format)
  set.seed(seed)
  n     <- nrow(geno_matrix)
  p     <- ncol(geno_matrix)
  qtl   <- sample(p, min(n_qtl, p))
  effs  <- rnorm(length(qtl), 0, 1)
  g_val <- as.numeric(geno_matrix[, qtl, drop = FALSE] %*% effs)
  # Scale genetic variance to h2, residual to (1-h2)
  var_g <- var(g_val)
  if (var_g < 1e-10) var_g <- 1
  e_val <- rnorm(n, 0, sqrt(var_g * (1 - h2) / h2))
  blues <- round((g_val + e_val - mean(g_val + e_val)) / sd(g_val + e_val), 4)
  names(blues) <- rownames(geno_matrix)
  if (identical(format, "data.frame")) {
    return(data.frame(id = names(blues), blue = unname(blues),
                      stringsAsFactors = FALSE))
  }
  blues
}

# -- Flat-file writers ---------------------------------------------------------

#' Simulate a random phased haplotype pair
#' @param n Number of individuals. Default 40.
#' @param p Number of SNPs. Default 20.
#' @param seed RNG seed for reproducibility.
#' @return List with hap1, hap2, dosage (SNPs x individuals, values 0/1 for
#'   hap1/hap2, 0/1/2 for dosage = hap1 + hap2; numeric matrices), and
#'   sample_ids (character vector). Row names = SNP IDs (matching
#'   make_snpinfo(p)$SNP), column names = individual IDs. Matches
#'   read_phased_vcf()'s output shape exactly (see R/haplotypes.R, its
#'   final `list(hap1=, hap2=, dosage=, snp_info=, sample_ids=, phased=,
#'   phase_method=)`), so it's a drop-in for functions documented as taking
#'   read_phased_vcf() output (extract_haplotypes(), ga_vs_ts_simulation(),
#'   etc.). extract_haplotypes()'s phased-list branch specifically requires
#'   `dosage` and `sample_ids` (not just hap1/hap2) to determine dosage-
#'   matrix orientation -- omitting them makes it fall through to a
#'   NULL-`dosage` code path that errors with "argument is of length zero"
#'   (a real mismatch caught by running the test suite, not by inspection).
make_phased <- function(n = 40, p = 20, seed = 11L) {
  set.seed(seed)
  h1 <- matrix(sample(0:1, n * p, replace = TRUE), nrow = p, ncol = n)
  h2 <- matrix(sample(0:1, n * p, replace = TRUE), nrow = p, ncol = n)
  rownames(h1) <- rownames(h2) <- paste0("rs", seq_len(p))
  colnames(h1) <- colnames(h2) <- paste0("ind", seq_len(n))
  ids <- colnames(h1)
  list(hap1 = h1, hap2 = h2, dosage = h1 + h2, sample_ids = ids)
}

#' Write a numeric dosage CSV compatible with read_geno(format = "numeric")
write_numeric_csv <- function(G, info, sep = ",") {
  if (!"REF" %in% names(info)) info$REF <- "A"
  if (!"ALT" %in% names(info)) info$ALT <- "T"
  meta <- info[, c("SNP","CHR","POS","REF","ALT")]
  df   <- cbind(meta, as.data.frame(t(G)))
  path <- tempfile(fileext = ".csv")
  write.table(df, path, sep = sep, row.names = FALSE, quote = FALSE)
  path
}

#' Write a VCF file compatible with read_geno(format = "vcf")
write_vcf <- function(G, info) {
  gt_enc <- function(g) vapply(g, function(x) {
    if (is.na(x)) return("./.")
    c("0" = "0/0", "1" = "0/1", "2" = "1/1")[as.character(x)]
  }, character(1))
  path  <- tempfile(fileext = ".vcf")
  lines <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",  # required by SNPRelate
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
            rownames(G)), collapse = "\t")
  )
  ref <- if ("REF" %in% names(info)) info$REF else rep("A", nrow(info))
  alt <- if ("ALT" %in% names(info)) info$ALT else rep("T", nrow(info))
  for (i in seq_len(nrow(info))) {
    lines <- c(lines, paste(c(
      info$CHR[i], info$POS[i], info$SNP[i],
      ref[i], alt[i], ".", "PASS", ".", "GT",
      gt_enc(G[, i])
    ), collapse = "\t"))
  }
  writeLines(lines, path)
  path
}

#' Write a HapMap file compatible with read_geno(format = "hapmap")
write_hapmap <- function(G, info) {
  ref <- if ("REF" %in% names(info)) info$REF else rep("A", nrow(info))
  alt <- if ("ALT" %in% names(info)) info$ALT else rep("T", nrow(info))
  hmp_decode <- function(g, r, a) {
    vapply(g, function(x) {
      if (is.na(x)) return("NN")
      c("0" = paste0(r, r), "1" = paste0(r, a),
        "2" = paste0(a, a))[as.character(x)]
    }, character(1))
  }
  sample_calls <- mapply(
    function(r, a, col) hmp_decode(col, r, a),
    ref, alt, as.list(as.data.frame(t(G))),
    SIMPLIFY = FALSE
  )
  sample_mat <- t(do.call(cbind, sample_calls))
  hdr <- data.frame(
    "rs#"       = info$SNP,
    alleles     = paste0(ref, "/", alt),
    chrom       = info$CHR,
    pos         = info$POS,
    strand      = "+",
    "assembly#" = "NA",
    center      = "NA",
    protLSID    = "NA",
    assayLSID   = "NA",
    panelLSID   = "NA",
    QCcode      = "NA",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(hdr, sample_mat)
  colnames(out)[12:ncol(out)] <- rownames(G)
  path <- tempfile(fileext = ".hmp.txt")
  write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}
