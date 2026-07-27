# -----------------------------------------------------------------------------
# Big_LD.R  -  Core per-chromosome LD block segmentation
#              Supports standard r^2 (default) and kinship-adjusted rV^2
#              C++ kernels used throughout for speed
# -----------------------------------------------------------------------------

#' @keywords internal
# -----------------------------------------------------------------------------
# LD-informed overlap resolution
# -----------------------------------------------------------------------------
#' @noRd
.resolve_overlap <- function(blocks,
                             adj_mat,
                             k_rep = 10L,
                             singleton_as_block = FALSE) {

  if (is.null(blocks) || nrow(blocks) == 0L) {
    return(blocks)
  }

  blocks <- as.matrix(blocks)
  storage.mode(blocks) <- "integer"

  # Keep only valid intervals. Invalid intervals can appear after boundary
  # correction if a block is completely absorbed by its neighbour.
  blocks <- blocks[
    is.finite(blocks[, 1L]) &
      is.finite(blocks[, 2L]) &
      blocks[, 1L] <= blocks[, 2L],
    ,
    drop = FALSE
  ]

  if (!nrow(blocks)) {
    return(blocks)
  }

  # Deterministic ordering is essential. The overlap algorithm is local:
  # it compares block i with block i+1, so intervals must be sorted by
  # start position and then by end position.
  blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

  # If singleton_as_block = FALSE, singleton intervals are not valid output
  # blocks. They are removed before overlap resolution because they do not
  # provide a left/right LD core and can prevent stable interval partitioning.
  #
  # If singleton_as_block = TRUE, singleton intervals are retained at this
  # stage, but if a singleton is fully contained inside a larger block it is
  # removed because the marker is already represented by the multi-SNP block.
  singleton_idx <- which(blocks[, 1L] == blocks[, 2L])

  if (length(singleton_idx)) {

    remove_singleton <- integer(0L)

    for (si in singleton_idx) {

      pos <- blocks[si, 1L]

      host <- which(
        seq_len(nrow(blocks)) != si &
          blocks[, 1L] <= pos &
          blocks[, 2L] >= pos
      )

      if (!isTRUE(singleton_as_block) || length(host) > 0L) {
        remove_singleton <- c(remove_singleton, si)
      }
    }

    if (length(remove_singleton)) {
      blocks <- blocks[-unique(remove_singleton), , drop = FALSE]
    }
  }

  if (!nrow(blocks)) {
    return(blocks)
  }

  blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

  # Remove fully contained non-singleton intervals before pairwise splitting.
  # Example:
  #   A = [20, 100]
  #   B = [40, 60]
  # Here B carries no unique SNP-index interval outside A. There is no
  # meaningful left/right boundary to estimate; B is redundant and should not
  # survive as a separate LD block.
  repeat {

    n0 <- nrow(blocks)
    keep <- rep(TRUE, n0)

    for (i in seq_len(n0)) {

      if (!keep[i]) next

      contained <- which(
        keep &
          seq_len(n0) != i &
          blocks[, 1L] >= blocks[i, 1L] &
          blocks[, 2L] <= blocks[i, 2L]
      )

      if (length(contained)) {
        keep[contained] <- FALSE
      }
    }

    blocks <- blocks[keep, , drop = FALSE]
    blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

    if (nrow(blocks) == n0) break
  }

  if (nrow(blocks) < 2L) {
    return(blocks)
  }

  i <- 1L

  while (i < nrow(blocks)) {

    sA <- blocks[i, 1L]
    eA <- blocks[i, 2L]
    sB <- blocks[i + 1L, 1L]
    eB <- blocks[i + 1L, 2L]

    if (eA < sB) {
      i <- i + 1L
      next
    }

    # Three zones, defined by SNP index:
    #
    #   Block A: [sA ........ sB-1 | sB ........ eA]
    #   Block B:              [sB ........ eA | eA+1 ........ eB]
    #
    #   left_core  = A-exclusive SNPs
    #   overlap    = disputed SNPs present in both A and B
    #   right_core = B-exclusive SNPs
    #
    # seq.int() is used with explicit conditions to avoid R's descending
    # sequence behaviour when start > end.
    left_core <- if (sB > sA) {
      seq.int(sA, sB - 1L)
    } else {
      integer(0L)
    }

    overlap <- if (eA >= sB) {
      seq.int(sB, eA)
    } else {
      integer(0L)
    }

    right_core <- if (eB > eA) {
      seq.int(eA + 1L, eB)
    } else {
      integer(0L)
    }

    # Fallback: if one side has no exclusive core, one interval is fully
    # contained in the other or both start at the same position. There is no
    # valid left/right LD anchor for a boundary decision.
    #
    # In that case, retain the union as one block and remove the redundant
    # neighbour. This is safer than preserving nested intervals because nested
    # intervals cause duplicate GWAS/QTL assignments.
    has_left <- length(left_core) > 0L
    has_right <- length(right_core) > 0L

    if (!has_left || !has_right) {

      blocks[i, ] <- c(min(sA, sB), max(eA, eB))
      blocks <- blocks[-(i + 1L), , drop = FALSE]
      blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

      # Recheck the previous position because the union may now overlap
      # the block before it.
      i <- max(1L, i - 1L)
      next
    }

    # Representatives are sampled from the boundary-adjacent core regions.
    #
    # Left representatives:
    #   last k_rep SNPs of left_core, i.e. closest SNPs before the overlap.
    #
    # Right representatives:
    #   first k_rep SNPs of right_core, i.e. closest SNPs after the overlap.
    #
    # This keeps the boundary decision local and prevents distant SNPs in
    # large blocks from dominating the assignment.
    k_L <- min(k_rep, length(left_core))
    k_R <- min(k_rep, length(right_core))

    left_reps <- tail(left_core, k_L)
    right_reps <- head(right_core, k_R)

    # Compute signed LD evidence for each disputed SNP:
    #
    #   score[s] = mean r2(s, left_reps) - mean r2(s, right_reps)
    #
    # Positive score:
    #   disputed SNP is more strongly associated with the left block.
    #
    # Negative score:
    #   disputed SNP is more strongly associated with the right block.
    #
    # NA values are treated as zero. This can occur for monomorphic or
    # numerically constant adjusted columns. A zero score is neutral and
    # does not force an artificial boundary.
    scores <- vapply(
      overlap,
      function(s) {
        r2v <- as.numeric(col_r2_cpp(adj_mat, s))

        r2_L <- mean(r2v[left_reps], na.rm = TRUE)
        r2_R <- mean(r2v[right_reps], na.rm = TRUE)

        if (is.na(r2_L) || is.na(r2_R)) {
          0
        } else {
          r2_L - r2_R
        }
      },
      numeric(1L)
    )

    # Cumulative-score boundary rule.
    #
    # Instead of assigning each disputed SNP independently, we cumulate the
    # signed evidence:
    #
    #   cum_score[k] = sum(scores[1:k])
    #
    # The split is placed at the last SNP where cumulative support still
    # favours the left block. This prevents local noisy scores from creating
    # jagged or biologically implausible block boundaries.
    #
    # Example:
    #   scores    = [+0.3, -0.1, +0.5, -0.3, -0.5, -0.7]
    #   cum_score = [+0.3, +0.2, +0.7, +0.4, -0.1, -0.8]
    #
    # The last non-negative cumulative score is position 4, so SNPs 1-4
    # are assigned to the left block and SNPs 5-6 to the right block.
    #
    # All-negative case:
    #   scores    = [-0.3, -0.5, -0.2]
    #   cum_score = [-0.3, -0.8, -1.0]
    #
    # No position has cum_score >= 0, so all disputed SNPs go to the right
    # block.
    #
    # All-zero case:
    #   scores    = [0, 0, 0]
    #   cum_score = [0, 0, 0]
    #
    # All positions satisfy cum_score >= 0. The overlap is assigned to the
    # left block, consistent with keeping the first-detected block in a
    # complete tie.
    cum_score <- cumsum(scores)
    last_left <- max(c(0L, which(cum_score >= 0)))

    if (last_left == 0L) {

      # Every disputed SNP belongs to block B.
      # Block A shrinks to [sA .. sB-1].
      # Block B keeps [sB .. eB].
      blocks[i, ] <- c(sA, sB - 1L)
      blocks[i + 1L, ] <- c(sB, eB)

    } else if (last_left == length(overlap)) {

      # Every disputed SNP belongs to block A.
      # Block A keeps [sA .. eA].
      # Block B shrinks to [eA+1 .. eB].
      #
      # This may create an invalid or singleton block B if B had little
      # or no right-side support. The clean-up step below removes such
      # intervals when singleton_as_block = FALSE.
      blocks[i, ] <- c(sA, eA)
      blocks[i + 1L, ] <- c(eA + 1L, eB)

    } else {

      # Mixed assignment:
      # SNPs up to split_snp go to A; subsequent SNPs go to B.
      split_snp <- overlap[last_left]

      blocks[i, ] <- c(sA, split_snp)
      blocks[i + 1L, ] <- c(split_snp + 1L, eB)
    }

    # General clean-up after the split:
    #   1. Remove invalid intervals with start > end.
    #   2. Remove singletons when singleton_as_block = FALSE.
    #   3. Remove contained intervals again, because a split/union can create
    #      a redundant interval relative to its neighbour.
    blocks <- blocks[
      blocks[, 1L] <= blocks[, 2L],
      ,
      drop = FALSE
    ]

    if (!isTRUE(singleton_as_block)) {
      blocks <- blocks[
        blocks[, 1L] < blocks[, 2L],
        ,
        drop = FALSE
      ]
    }

    if (nrow(blocks) > 1L) {

      blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

      keep <- rep(TRUE, nrow(blocks))

      for (ii in seq_len(nrow(blocks))) {

        if (!keep[ii]) next

        contained <- which(
          keep &
            seq_len(nrow(blocks)) != ii &
            blocks[, 1L] >= blocks[ii, 1L] &
            blocks[, 2L] <= blocks[ii, 2L]
        )

        if (length(contained)) {
          keep[contained] <- FALSE
        }
      }

      blocks <- blocks[keep, , drop = FALSE]
    }

    blocks <- blocks[order(blocks[, 1L], blocks[, 2L]), , drop = FALSE]

    if (nrow(blocks) < 2L) {
      break
    }

    # Move one step back because the current boundary update may have created
    # a new overlap with the preceding interval.
    i <- max(1L, i - 1L)
  }

  blocks
}

#' @keywords internal
# -----------------------------------------------------------------------------
# LD-informed gap resolution
# -----------------------------------------------------------------------------
#' @noRd
.close_interblock_gaps <- function(blocks, adj_mat, k_rep = 10L) {

  if (is.null(blocks) || nrow(blocks) < 2L) {
    return(blocks)
  }

  blocks <- blocks[order(blocks[, 1L]), , drop = FALSE]

  for (i in seq_len(nrow(blocks) - 1L)) {

    sA <- blocks[i, 1L]
    eA <- blocks[i, 2L]
    sB <- blocks[i + 1L, 1L]
    eB <- blocks[i + 1L, 2L]

    # No SNP-index gap
    if ((eA + 1L) >= sB) {
      next
    }

    gap_snps <- seq.int(eA + 1L, sB - 1L)

    # Representatives close to the gap
    left_core <- seq.int(sA, eA)
    right_core <- seq.int(sB, eB)

    k_L <- min(k_rep, length(left_core))
    k_R <- min(k_rep, length(right_core))

    left_reps <- tail(left_core, k_L)
    right_reps <- head(right_core, k_R)

    scores <- vapply(
      gap_snps,
      function(s) {
        r2v <- as.numeric(col_r2_cpp(adj_mat, s))

        r2_L <- mean(r2v[left_reps], na.rm = TRUE)
        r2_R <- mean(r2v[right_reps], na.rm = TRUE)

        if (is.na(r2_L) || is.na(r2_R)) {
          return(0)
        }

        r2_L - r2_R
      },
      numeric(1L)
    )

    cum_score <- cumsum(scores)
    last_left <- max(c(0L, which(cum_score >= 0)))

    if (last_left == 0L) {

      # All gap SNPs assigned to right block
      blocks[i, 2L] <- eA
      blocks[i + 1L, 1L] <- eA + 1L

    } else if (last_left == length(gap_snps)) {

      # All gap SNPs assigned to left block
      blocks[i, 2L] <- sB - 1L
      blocks[i + 1L, 1L] <- sB

    } else {

      # Split gap between left and right block
      split_snp <- gap_snps[last_left]

      blocks[i, 2L] <- split_snp
      blocks[i + 1L, 1L] <- split_snp + 1L
    }
  }

  blocks
}

#' @keywords internal
# -----------------------------------------------------------------------------
# LD-informed absorption of unassigned singleton SNPs
# -----------------------------------------------------------------------------
#' @noRd
.absorb_unassigned_snps <- function(blocks,
                                    unassigned_idx,
                                    adj_mat,
                                    k_rep = 2L) {

  if (is.null(blocks) || nrow(blocks) == 0L || !length(unassigned_idx)) {
    return(blocks)
  }

  blocks <- as.matrix(blocks)
  storage.mode(blocks) <- "integer"

  blocks <- blocks[
    blocks[, 1L] <= blocks[, 2L],
    ,
    drop = FALSE
  ]

  if (!nrow(blocks)) {
    return(blocks)
  }

  blocks <- blocks[
    order(blocks[, 1L], blocks[, 2L]),
    ,
    drop = FALSE
  ]

  unassigned_idx <- sort(unique(as.integer(unassigned_idx)))

  for (s in unassigned_idx) {

    # Skip SNPs already covered by an LD block.
    if (any(blocks[, 1L] <= s & blocks[, 2L] >= s)) {
      next
    }

    left_candidates <- which(blocks[, 2L] < s)
    right_candidates <- which(blocks[, 1L] > s)

    left_i <- if (length(left_candidates)) {
      left_candidates[which.max(blocks[left_candidates, 2L])]
    } else {
      NA_integer_
    }

    right_i <- if (length(right_candidates)) {
      right_candidates[which.min(blocks[right_candidates, 1L])]
    } else {
      NA_integer_
    }

    if (!is.na(left_i) && !is.na(right_i)) {

      left_core <- seq.int(blocks[left_i, 1L], blocks[left_i, 2L])
      right_core <- seq.int(blocks[right_i, 1L], blocks[right_i, 2L])

      k_L <- min(k_rep, length(left_core))
      k_R <- min(k_rep, length(right_core))

      left_reps <- tail(left_core, k_L)
      right_reps <- head(right_core, k_R)

      r2v <- as.numeric(col_r2_cpp(adj_mat, s))

      r2_L <- mean(r2v[left_reps], na.rm = TRUE)
      r2_R <- mean(r2v[right_reps], na.rm = TRUE)

      if (is.na(r2_L)) r2_L <- 0
      if (is.na(r2_R)) r2_R <- 0

      if (r2_L >= r2_R) {
        blocks[left_i, 2L] <- max(blocks[left_i, 2L], s)
      } else {
        blocks[right_i, 1L] <- min(blocks[right_i, 1L], s)
      }

    } else if (!is.na(left_i)) {

      blocks[left_i, 2L] <- max(blocks[left_i, 2L], s)

    } else if (!is.na(right_i)) {

      blocks[right_i, 1L] <- min(blocks[right_i, 1L], s)
    }
  }

  blocks <- blocks[
    order(blocks[, 1L], blocks[, 2L]),
    ,
    drop = FALSE
  ]

  blocks
}

#' LD Block Segmentation (r^2 or rV^2, C++ accelerated)
#'
#' @description
#' Core per-chromosome LD block detection. Two LD metrics are supported:
#'
#' \describe{
#'   \item{\code{method = "r2"} (default)}{Standard squared Pearson
#'     correlation, computed by the C++ Armadillo kernel. No kinship matrix
#'     is required. Suitable for large datasets (hundreds of thousands to
#'     millions of markers) and unstructured or mildly structured populations.}
#'   \item{\code{method = "rV2"}}{Kinship-adjusted squared correlation
#'     (Kim et al. 2018). Requires computing and inverting a GRM via
#'     AGHmatrix + ASRgenomics. Recommended for highly related populations
#'     (livestock, inbred lines) with moderate marker counts (< 200 k per chr).}
#' }
#'
#' By default, SNP-index gaps between adjacent LD blocks are resolved using
#' an LD-informed boundary assignment procedure. This ensures that every SNP
#' passing quality-control filters is assigned to exactly one LD block,
#' preventing unassigned markers between neighbouring blocks while preserving
#' biologically meaningful physical gaps where no SNPs are present.
#'
#' For genome-wide analyses use \code{\link{run_Big_LD_all_chr}}.
#' For automatic parameter tuning use \code{\link{tune_LD_params}}.
#'
#' @param geno Numeric matrix (individuals x SNPs), 0/1/2. Row names used
#'   as individual IDs; auto-generated if absent.
#' @param SNPinfo Data frame: col 1 = rsID, col 2 = bp position.
#' @param method Character. \code{"r2"} (default) or \code{"rV2"}.
#' @param CLQcut Numeric in (0,1]. LD threshold for clique edges. Default 0.5.
#' @param clstgap Integer. Max bp gap within clique (split=TRUE). Default 40000.
#' @param leng Integer. Boundary-scan half-window (SNPs). Default 200.
#' @param subSegmSize Integer. Max SNPs per CLQD call. Default 1500.
#' @param MAFcut Numeric. Minor allele frequency minimum. Default 0.05.
#' @param appendrare Logical. Append rare SNPs after block detection. Default FALSE.
#' @param singleton_as_block Logical. If \code{TRUE}, every SNP that passes MAF
#'   filtering but has pairwise r^2 below \code{CLQcut} with all neighbours
#'   (i.e. is not assigned to any clique) is returned as a single-SNP block with
#'   \code{start == end} and \code{length_bp == 1}. Default \code{FALSE}.
#'   These blocks are excluded from haplotype analysis by the default
#'   \code{min_snps = 3} threshold in \code{extract_haplotypes()}, but are
#'   useful for auditing coverage and for single-SNP feature engineering.
#'
#' @param close_gaps_with_snps Logical. If \code{TRUE}, adjacent LD blocks
#'   separated only by unassigned SNPs are expanded so that every SNP passing
#'   MAF filtering is assigned to a block. Gap SNPs are allocated to the
#'   flanking LD blocks using an LD-informed boundary rule based on their
#'   relative correlation with representative SNPs from each neighbouring
#'   block. Physical gaps without SNPs remain unchanged. Default
#'   \code{TRUE}.
#'
#' @param gap_k_rep Integer. Number of representative SNPs sampled from each
#'   side of an inter-block gap when \code{close_gaps_with_snps = TRUE}.
#'   Representatives are chosen closest to the gap boundary and are used to
#'   evaluate whether each gap SNP is more strongly associated with the left
#'   or right neighbouring block. Larger values provide more stable boundary
#'   assignment at increased computational cost. Default \code{2L}.
#'
#' @param checkLargest Logical. Dense-core pre-pass for large windows. Default FALSE.
#' @param CLQmode \code{"Density"} (default) or \code{"Maximal"}.
#' @param kin_method Character. GRM whitening: \code{"chol"} (default) or
#'   \code{"eigen"}. Ignored when \code{method = "r2"}.
#' @param split Logical. Split cliques at gaps > clstgap. Default FALSE.
#' @param max_bp_distance Max bp for r\eqn{^2} (\code{0L} = all). Default \code{0L}
#' @param digits Integer. LD rounding: \code{-1} (default, no rounding)
#'   or a positive integer.
#' @param n_threads Integer. OpenMP threads for C++ LD kernel. Default 1.
#' @param seed Integer or NULL. Set for reproducibility.
#' @param verbose Logical. Progress messages.
#'
#' @return data.frame with columns: start, end, start.rsID, end.rsID,
#'   start.bp, end.bp.
#'
#' @seealso \code{\link{run_Big_LD_all_chr}}, \code{\link{CLQD}},
#'   \code{\link{compute_r2}}, \code{\link{compute_rV2}},
#'   \code{\link{tune_LD_params}}
#'
#' @references
#' Kim S-A et al. (2018) GENETICS 209(3):855-868.\cr
#' VanRaden PM (2008) J. Dairy Sci. 91(11):4414-4423.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' m <- 80; b1 <- 60; b2 <- 60
#' make_block <- function(m, size, p, flip = 0.03) {
#'   seed_snp <- rbinom(m, 2, p)
#'   M <- matrix(seed_snp, nrow = m, ncol = size)
#'   for (j in seq_len(size)) {
#'     idx <- sample.int(m, max(1L, floor(flip * m)))
#'     M[idx, j] <- pmin(2L, pmax(0L, M[idx, j] + sample(c(-1L, 1L), length(idx), TRUE)))
#'   }
#'   M
#' }
#' G <- cbind(make_block(m, b1, 0.30), make_block(m, b2, 0.60))
#' colnames(G) <- paste0("rs", seq_len(ncol(G)))
#' rownames(G) <- paste0("ind", seq_len(nrow(G)))
#' pos <- c(seq(1, by = 1000, length.out = b1), seq(5e6, by = 1000, length.out = b2))
#' SNPinfo <- data.frame(SNP = colnames(G), POS = pos)
#' blocks <- HapBlockR:::Big_LD(G, SNPinfo, method = "r2", CLQcut = 0.6,
#'                  leng = 30, subSegmSize = 120, verbose = FALSE)
#' head(blocks)
#' }
#' @export

Big_LD <- function(
    geno,
    SNPinfo,
    method       = c("r2", "rV2"),
    CLQcut       = 0.5,
    clstgap      = 40000,
    leng         = 200,
    subSegmSize  = 1500,
    MAFcut       = 0.05,
    appendrare   = FALSE,
    singleton_as_block = FALSE,
    close_gaps_with_snps = TRUE,
    gap_k_rep = 2L,
    checkLargest = FALSE,
    CLQmode         = c("Density", "Maximal", "Louvain", "Leiden"),
    kin_method      = "chol",
    split           = FALSE,
    max_bp_distance = 0L,
    digits       = -1L,
    n_threads    = 1L,
    seed         = NULL,
    verbose      = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)

  if (!is.matrix(geno)) geno <- as.matrix(geno)
  if (ncol(geno) != nrow(SNPinfo))
    stop("ncol(geno) [", ncol(geno), "] != nrow(SNPinfo) [", nrow(SNPinfo), "]")
  if (ncol(SNPinfo) != 2L)
    stop("SNPinfo must have exactly 2 columns: [rsID, bp_position]")

  ids <- rownames(geno)
  if (is.null(ids)) { ids <- sprintf("ind%04d", seq_len(nrow(geno))); rownames(geno) <- ids }

  # -- Prepare (centre or whiten) --------------------------------------------
  prep       <- prepare_geno(geno, method = method, kin_method = kin_method, verbose = verbose)
  adj_geno   <- prep$adj_geno
  V_inv_sqrt <- prep$V_inv_sqrt   # NULL for r^2

  Ogeno    <- geno
  OSNPinfo <- SNPinfo

  # -- MAF + monomorphic filter (C++) ---------------------------------------
  keep_poly <- maf_filter_cpp(geno, maf_cut = MAFcut)
  monoSNPs  <- OSNPinfo[!keep_poly, , drop = FALSE]
  geno      <- geno[, keep_poly, drop = FALSE]
  adjN      <- adj_geno[, keep_poly, drop = FALSE]
  SNPinfo   <- SNPinfo[keep_poly, , drop = FALSE]

  if (ncol(geno) < 2L) {
    warning("[Big_LD] Fewer than 2 polymorphic SNPs after MAF filter. Returning empty.")
    return(data.frame(start = integer(), end = integer(),
                      start.rsID = character(), end.rsID = character(),
                      start.bp = numeric(), end.bp = numeric(),
                      n_snps = integer()))
  }

  # ==========================================================================
  # Internal helpers
  # ==========================================================================

  # -- cutsequence.modi: find weak-LD boundaries using C++ boundary_scan ----
  cutsequence.modi <- function(geno, adj_geno, leng, subSegmSize, CLQcut, digits, n_threads) {
    p <- ncol(geno)
    if (p <= subSegmSize) return(list(p, NULL))
    if (isTRUE(verbose)) message("[Big_LD] Subsegmenting via C++ boundary scan...")

    modeNum <- 1L; lastnum <- 0L; cutpoints <- NULL; i <- leng

    while (i <= (p - leng)) {
      if ((i - lastnum) > 5L * subSegmSize) { modeNum <- 2L; break }

      found_cut <- FALSE
      for (j in unique(c(1L, min(10L, leng), leng))) {
        l_start <- i - j + 1L; l_end <- i
        r_start <- i + 1L;     r_end <- i + j
        if (r_end > p) next

        # C++ cross-boundary r^2 check
        sub_cols  <- l_start:r_end
        sub_adj   <- adj_geno[, sub_cols, drop = FALSE]
        cross_r2  <- compute_r2_cpp(sub_adj,
                                    digits    = as.integer(digits),
                                    n_threads = as.integer(n_threads))
        nl <- j; nr <- j
        cross_block <- cross_r2[seq_len(nl), (nl + 1L):(nl + nr), drop = FALSE]
        cross_block[cross_block < CLQcut] <- 0

        if (sum(cross_block, na.rm = TRUE) > 0) { i <- i + 1L; break }
        if (j == leng) found_cut <- TRUE
      }

      if (found_cut) {
        cutpoints <- c(cutpoints, i)
        lastnum   <- i
        if (isTRUE(verbose)) message(sprintf("[Big_LD]   Cut at SNP %d", i))
        i <- i + floor(leng / 2L)
      }
    }

    if (modeNum == 1L) {
      cutpoints <- c(0L, cutpoints, p)
      atfcut    <- NULL
      while (max(diff(cutpoints)) > subSegmSize) {
        diffseq    <- diff(cutpoints)
        recutpoint <- which(diffseq > subSegmSize)
        tt         <- cbind(cutpoints[recutpoint] + 1L, cutpoints[recutpoint + 1L])
        numvec     <- NULL
        tick       <- max(1L, as.integer(leng / 5L))

        for (k in seq_len(nrow(tt))) {
          st <- tt[k, 1L]; ed <- tt[k, 2L]
          if (ed > (p - leng)) ed <- p - leng
          cands <- (st + leng):(ed - leng)
          if (length(cands) == 0L) next

          # C++ boundary scan over candidate positions
          scan_cols <- st:(ed + leng)
          scan_cols <- scan_cols[scan_cols >= 1L & scan_cols <= p]
          scan_adj  <- adj_geno[, scan_cols, drop = FALSE]
          weak_flags <- boundary_scan_cpp(scan_adj,
                                          start     = tick + 1L,
                                          end       = ncol(scan_adj) - tick,
                                          half_w    = tick,
                                          threshold = CLQcut)
          # find positions with fewest cross-LD (weak_flags == 1)
          weakcount    <- as.integer(!weak_flags)   # 0 = weak, good for cut
          weakcount.s  <- sort(weakcount)
          weaks        <- weakcount.s[min(10L, length(weakcount.s))]
          rel_idx      <- which(weakcount <= weaks)
          abs_idx      <- rel_idx + st + leng - 1L
          abs_idx      <- abs_idx[abs_idx >= st & abs_idx <= ed]
          if (length(abs_idx) == 0L) abs_idx <- as.integer(round((st + ed) / 2L))
          nearcenter   <- vapply(abs_idx, function(x) abs((ed - x) - (x - st)), numeric(1L))
          addcut       <- abs_idx[which.min(nearcenter)]
          numvec       <- c(numvec, addcut)
          atfcut       <- c(atfcut, addcut)
        }
        cutpoints <- sort(c(cutpoints, numvec))
        if (!any(diff(cutpoints) > subSegmSize)) break
      }
    } else {
      cutpoints <- seq(subSegmSize, p, subSegmSize / 2L)
      atfcut    <- cutpoints
      if (max(cutpoints) == p) atfcut <- atfcut[-length(atfcut)]
      else                     cutpoints <- c(cutpoints, p)
    }
    list(cutpoints, atfcut)
  }

  # -- intervalCliqueList ----------------------------------------------------
  intervalCliqueList <- function(clstlist, allsnps, onlybp) {
    bp.clstlist <- lapply(clstlist, function(x) onlybp[x])
    bp.allsnps  <- lapply(allsnps,  function(x) onlybp[x])
    IMsize      <- length(bp.clstlist)
    adjacencyM  <- matrix(0L, IMsize, IMsize)
    for (i in seq_len(IMsize))
      for (j in seq_len(IMsize))
        adjacencyM[i, j] <- length(intersect(bp.allsnps[[i]], bp.allsnps[[j]]))
    diag(adjacencyM) <- 0L
    ig <- igraph::graph_from_adjacency_matrix(adjacencyM, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
    icliques <- if (max(igraph::coreness(ig)) > 10L) igraph::max_cliques(ig, min = 1L)
    else                                  igraph::cliques(ig, min = 1L)
    icliques  <- icliques[order(vapply(icliques, min, numeric(1L)))]
    intervals  <- lapply(icliques, function(x) sort(unique(unlist(bp.clstlist[x]))))
    weight.itv <- vapply(intervals, length, integer(1L))
    ir         <- t(vapply(intervals, range, numeric(2L)))
    ur         <- unique(ir)
    ri         <- cbind(ir, weight.itv)
    info_idx   <- apply(ur, 1L, function(x) {
      same <- intersect(which(ri[,1L] == x[1L]), which(ri[,2L] == x[2L]))
      same[which.max(ri[same, 3L])]
    })
    list(intervals[info_idx], ri[info_idx, 3L])
  }

  # -- find.maximum.indept ---------------------------------------------------
  find.maximum.indept <- function(sample.itv, sample.weight) {
    n              <- length(sample.itv)
    interval.range <- t(vapply(sample.itv, range, numeric(2L)))
    pre.range      <- vector("list", n)
    for (i in seq_len(n)) {
      idx <- which(interval.range[, 2L] < interval.range[i, 1L])
      pre.range[[i]] <- if (length(idx) > 0L) idx else NA_integer_
    }
    sources <- which(vapply(pre.range, function(x) all(is.na(x)), logical(1L)))
    if (length(sources) < n) {
      not.s <- setdiff(seq_len(n), sources)
      for (i in not.s) {
        pp <- sort(unique(unlist(pre.range[pre.range[[i]]])))
        pre.range[[i]] <- setdiff(pre.range[[i]], pp)
      }
      rw  <- rep(0, n); rw[sources] <- sample.weight[sources]
      ptr <- rep(0L, n); ptr[sources] <- NA_integer_
      info <- cbind(seq_len(n), rw, ptr, expl = replace(rep(0L,n), sources, 1L))
      for (i in not.s) {
        mp   <- pre.range[[i]]
        ni   <- info[mp, , drop = FALSE]
        mi   <- ni[ni[,2L] == max(ni[,2L]), , drop = FALSE][1L,]
        info[i,2L] <- sample.weight[i] + mi[2L]
        info[i,3L] <- mi[1L]; info[i,4L] <- 1L
      }
      si  <- which(info[,2L] == max(info[,2L]))[1L]
      pre <- info[si,3L]; is <- c(pre, si)
      while (!is.na(pre)) { pre <- info[pre,3L]; is <- c(pre, is) }
      list(indept.set = is[!is.na(is)], indept.set.weight = max(info[,2L]))
    } else {
      list(indept.set = which(sample.weight == max(sample.weight)),
           indept.set.weight = max(sample.weight))
    }
  }

  # -- constructLDblock ------------------------------------------------------
  constructLDblock <- function(clstlist, subSNPinfo) {
    safe_range <- function(v) {
      y <- as.numeric(stats::na.omit(v))
      if (length(y) < 2L) return(NULL)
      st <- min(y); ed <- max(y)
      if (!is.finite(st) || !is.finite(ed) || st > ed) return(NULL)
      c(start = st, end = ed)
    }
    SP   <- as.numeric(as.character(subSNPinfo[[2L]]))
    clst <- Filter(function(x) length(x) >= 2L,
                   lapply(clstlist, function(x) as.numeric(stats::na.omit(unlist(x)))))
    if (length(clst) == 0L) return(NULL)
    Total <- NULL
    while (length(clst) > 0L) {
      cr   <- Filter(Negate(is.null), lapply(clst, safe_range))
      asnp <- lapply(cr, function(rg) seq.int(rg["start"], rg["end"]))
      if (length(asnp) == 0L) break
      cand <- intervalCliqueList(clst, asnp, SP)
      ivs  <- cand[[1L]]; wts <- cand[[2L]]
      if (length(ivs) == 0L || length(wts) == 0L) break
      MWIS <- find.maximum.indept(ivs, wts)
      iset <- ivs[MWIS[[1L]]]
      sblk <- Filter(Negate(is.null), lapply(iset, function(x) {
        safe_range(as.numeric(stats::na.omit(match(x, SP))))
      }))
      if (length(sblk) == 0L) break
      sb   <- do.call(rbind, sblk)
      Total <- rbind(Total, sb)
      used  <- unique(unlist(lapply(sblk, function(rg) seq.int(rg["start"], rg["end"]))))
      clst  <- Filter(function(x) length(x) >= 2L, lapply(clst, function(x) setdiff(x, used)))
    }
    Total
  }

  # -- subBigLD --------------------------------------------------------------
  subBigLD <- function(sg, si, ag, CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads, max_bp_distance = 0L) {
    bv   <- CLQD(sg, si, ag, CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                 codechange = FALSE, checkLargest = checkLargest, split = split,
                 digits = digits, n_threads = n_threads, max_bp_distance = max_bp_distance, verbose = FALSE)
    if (all(is.na(bv))) return(matrix(integer(0), nrow=0L, ncol=2L))
    bins <- seq_len(max(bv[!is.na(bv)]))
    cl   <- lapply(bins, function(x) sort(which(bv == x)))
    cl   <- Filter(length, cl)  # drop empty bins (Louvain IDs may not be 1..n)
    if (length(cl)) cl <- cl[order(vapply(cl, min, numeric(1L)))]
    nowLD <- constructLDblock(cl, si)
    if (is.null(nowLD)) return(matrix(integer(0), nrow=0L, ncol=2L))
    if (!is.matrix(nowLD)) nowLD <- matrix(nowLD, nrow=1L)
    nowLD[order(nowLD[,1L]), , drop = FALSE]
  }

  # -- appendSGTs ------------------------------------------------------------
  appendSGTs <- function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap,
                         checkLargest, CLQmode, prep_full, split, digits, n_threads,
                         max_bp_distance = 0L) {
    if (isTRUE(verbose)) message("[Big_LD] appendSGTs: assigning rare SNPs.")
    expandB <- NULL

    # Recompute adjusted geno for full original set
    Ogeno_c <- scale(Ogeno, center = TRUE, scale = FALSE)
    ag_full <- if (!is.null(prep_full$V_inv_sqrt)) {
      prep_full$V_inv_sqrt %*% Ogeno_c
    } else {
      Ogeno_c
    }

    snp1 <- which(OSNPinfo[,2L] < LDblocks[1L, 5L])
    if (length(snp1) > 2L) {
      OSNPs      <- seq_len(max(snp1))
      firstB     <- LDblocks[1L, ]
      secondSNPs <- which(OSNPinfo[,2L] >= firstB$start.bp & OSNPinfo[,2L] <= firstB$end.bp)
      blk_adj    <- ag_full[, c(secondSNPs, OSNPs), drop = FALSE]
      rv2        <- compute_r2_cpp(blk_adj, digits = as.integer(digits), n_threads = as.integer(n_threads))
      rv2        <- rv2[seq_along(secondSNPs), (length(secondSNPs)+1L):ncol(rv2), drop = FALSE]
      cor2ratio  <- colSums(rv2 > CLQcut) / length(secondSNPs)
      cor2numT   <- c(cor2ratio > 0.6, 1L)
      points2    <- min(which(cor2numT > 0))
      NsecondSNPs <- points2:max(secondSNPs)
      reOSNPs    <- setdiff(seq_len(max(NsecondSNPs)), NsecondSNPs)
      if (length(reOSNPs) > 1L) {
        sb <- subBigLD(Ogeno[,reOSNPs], OSNPinfo[reOSNPs,], ag_full[,reOSNPs,drop=FALSE],
                       CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
        expandB <- rbind(expandB, sb + min(reOSNPs) - 1L)
      }
      firstSNPs <- NsecondSNPs
    } else {
      firstB    <- LDblocks[1L, ]
      firstSNPs <- which(OSNPinfo[,2L] >= firstB$start.bp & OSNPinfo[,2L] <= firstB$end.bp)
    }

    if (nrow(LDblocks) > 1L) {
      for (i in seq_len(nrow(LDblocks) - 1L)) {
        if (!length(firstSNPs)) next  # skip if first block has no SNPs in OSNPinfo
        secondB    <- LDblocks[i+1L, ]
        secondSNPs <- which(OSNPinfo[,2L] >= secondB$start.bp & OSNPinfo[,2L] <= secondB$end.bp)
        if (!length(secondSNPs)) {
          expandB <- rbind(expandB, range(firstSNPs))
          firstSNPs <- firstSNPs  # no next block found; keep and skip
          next
        }
        OSNPs      <- setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))
        if (length(OSNPs) == 0L) {
          expandB <- rbind(expandB, range(firstSNPs)); firstSNPs <- secondSNPs
        } else {
          rv1 <- compute_r2_cpp(ag_full[,c(firstSNPs,OSNPs),drop=FALSE],
                                digits=as.integer(digits), n_threads=as.integer(n_threads))
          rv1 <- rv1[seq_along(firstSNPs),(length(firstSNPs)+1L):ncol(rv1),drop=FALSE]
          rv2 <- compute_r2_cpp(ag_full[,c(secondSNPs,OSNPs),drop=FALSE],
                                digits=as.integer(digits), n_threads=as.integer(n_threads))
          rv2 <- rv2[seq_along(secondSNPs),(length(secondSNPs)+1L):ncol(rv2),drop=FALSE]
          cor1ratio <- colSums(rv1 > CLQcut) / length(firstSNPs)
          cor2ratio <- colSums(rv2 > CLQcut) / length(secondSNPs)
          cor1numT  <- c(1L, cor1ratio > 0.6, 0L)
          cor2numT  <- c(0L, cor2ratio > 0.6, 1L)
          p1  <- min(max(firstSNPs) + max(which(cor1numT > 0)) - 1L, ncol(Ogeno))
          p2  <- min(max(firstSNPs) + max(which(cor2numT > 0)) - 1L, ncol(Ogeno))
          p2  <- max(p2, min(secondSNPs))  # ensure p2 <= max(secondSNPs) direction
          NF  <- min(firstSNPs):p1
          NS  <- min(p2, max(secondSNPs)):max(secondSNPs)
          if (max(NF) < min(NS)) {
            expandB <- rbind(expandB, range(NF))
            reO <- setdiff(min(NF):max(NS), c(NF, NS))
            if (length(reO) > 1L) {
              sb <- subBigLD(Ogeno[,reO], OSNPinfo[reO,], ag_full[,reO,drop=FALSE],
                             CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
              if (nrow(sb) > 0L) expandB <- rbind(expandB, sb + min(reO) - 1L)
            }
            firstSNPs <- NS
          } else {
            rng <- min(firstSNPs):max(secondSNPs)
            sb  <- subBigLD(Ogeno[,rng], OSNPinfo[rng,], ag_full[,rng,drop=FALSE],
                            CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
            if (nrow(sb) == 0L) { firstSNPs <- secondSNPs; next }
            sb  <- sb + min(rng) - 1L
            if (nrow(sb) == 1L) firstSNPs <- sb[1L,1L]:sb[1L,2L]
            else { expandB <- rbind(expandB, sb[-nrow(sb),]); firstSNPs <- sb[nrow(sb),1L]:sb[nrow(sb),2L] }
          }
        }
      }
    }

    if (length(firstSNPs) > 0L && max(firstSNPs) < (ncol(Ogeno) - 1L)) {
      OSNPs <- (max(firstSNPs)+1L):ncol(Ogeno)
      rv2   <- compute_r2_cpp(ag_full[,c(firstSNPs,OSNPs),drop=FALSE],
                              digits=as.integer(digits), n_threads=as.integer(n_threads))
      rv2   <- rv2[seq_along(firstSNPs),(length(firstSNPs)+1L):ncol(rv2),drop=FALSE]
      cor1ratio <- colSums(rv2 > CLQcut) / length(firstSNPs)
      cor1numT  <- c(1L, cor1ratio > 0.6, 0L)
      p1 <- max(firstSNPs) + max(which(cor1numT > 0)) - 1L
      NF <- min(firstSNPs):p1
      expandB <- rbind(expandB, range(NF))
      reO <- setdiff(min(NF):ncol(Ogeno), NF)
      if (length(reO) > 1L) {
        sb <- subBigLD(Ogeno[,reO], OSNPinfo[reO,], ag_full[,reO,drop=FALSE],
                       CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
        expandB <- rbind(expandB, sb + min(reO) - 1L)
      }
    } else if (length(firstSNPs) > 0L) {
      expandB <- rbind(expandB, range(firstSNPs))
    }

    if (is.null(expandB) || nrow(expandB) == 0L) return(LDblocks)
    expandB <- expandB[expandB[,1L] != expandB[,2L], , drop = FALSE]
    data.frame(start      = expandB[,1L],
               end        = expandB[,2L],
               start.rsID = as.character(OSNPinfo[[1L]][expandB[,1L]]),
               end.rsID   = as.character(OSNPinfo[[1L]][expandB[,2L]]),
               start.bp   = as.numeric(OSNPinfo[[2L]][expandB[,1L]]),
               end.bp     = as.numeric(OSNPinfo[[2L]][expandB[,2L]]))
  }

  # ==========================================================================
  # Main loop
  # ==========================================================================

  cuts      <- cutsequence.modi(geno, adjN, leng, subSegmSize, CLQcut, digits, n_threads)
  cutpoints <- setdiff(cuts[[1L]], 0L)
  atfcut    <- cuts[[2L]]
  if (!is.null(atfcut)) atfcut <- sort(atfcut)

  cutblock <- cbind(c(1L, cutpoints + 1L), c(cutpoints, ncol(geno)))
  if (nrow(cutblock) > 1L) cutblock <- cutblock[-nrow(cutblock), , drop = FALSE]

  LDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2L)
  ld_count <- 0L  # running count of filled rows (avoids O(n) sum(!is.na()) per segment)

  for (i in seq_len(nrow(cutblock))) {

    nowst <- cutblock[i, 1L]
    nowed <- cutblock[i, 2L]

    sg <- geno[, nowst:nowed, drop = FALSE]
    ag <- adjN[, nowst:nowed, drop = FALSE]
    si <- SNPinfo[nowst:nowed, , drop = FALSE]

    bv <- CLQD(
      sg,
      si,
      ag,
      CLQcut = CLQcut,
      clstgap = clstgap,
      CLQmode = CLQmode,
      codechange = FALSE,
      checkLargest = checkLargest,
      split = split,
      digits = digits,
      n_threads = n_threads,
      max_bp_distance = max_bp_distance,
      verbose = verbose
    )

    # Collect unassigned SNP indices.
    #
    # NA in the CLQD bin vector means the SNP was not assigned to any LD clique.
    # These SNPs must always be tracked, irrespective of singleton_as_block.
    #
    # Later handling:
    #   singleton_as_block = TRUE
    #     -> keep genuinely isolated SNPs as single-SNP blocks.
    #
    #   singleton_as_block = FALSE
    #     -> absorb unassigned SNPs into the nearest flanking LD block.
    sing_local <- which(is.na(bv))
    sing_global <- sing_local + (nowst - 1L)

    if (!exists("singleton_idx_list")) {
      singleton_idx_list <- vector("list", nrow(cutblock))
    }

    singleton_idx_list[[i]] <- sing_global

    non_empty <- !all(is.na(bv))

    if (!non_empty) {

      nowLD <- matrix(
        integer(0),
        nrow = 0L,
        ncol = 2L
      )

    } else {

      bins <- seq_len(max(bv[!is.na(bv)]))

      cl <- lapply(
        bins,
        function(x) sort(which(bv == x))
      )

      # Drop empty bins. This is required because Louvain/Leiden/community
      # identifiers are not guaranteed to be contiguous integers.
      cl <- Filter(length, cl)

      if (length(cl)) {
        cl <- cl[order(vapply(cl, min, numeric(1L)))]
      }

      nowLD <- constructLDblock(cl, si)

      if (!is.null(nowLD)) {

        if (!is.matrix(nowLD)) {
          nowLD <- matrix(nowLD, nrow = 1L)
        }

        # Convert segment-local SNP indices to chromosome-local SNP indices.
        nowLD <- nowLD + (nowst - 1L)

        nowLD <- nowLD[
          order(nowLD[, 1L]),
          ,
          drop = FALSE
        ]

      } else {

        nowLD <- matrix(
          integer(0),
          nrow = 0L,
          ncol = 2L
        )
      }
    }

    if (nrow(nowLD) > 0L) {
      LDblocks[(ld_count + 1L):(ld_count + nrow(nowLD)), ] <- nowLD
      ld_count <- ld_count + nrow(nowLD)
    }

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[Big_LD] Segment %d/%d done.",
          i,
          nrow(cutblock)
        )
      )
    }
  }

  done <- LDblocks[seq_len(ld_count), , drop = FALSE]

  # -- Optional re-merge across forced cut-points ----------------------------
  # Skip re-merge when modeNum=2 (all cuts forced = no natural LD boundaries).
  # In modeNum=2, atfcut covers every cut position; re-running CLQD across
  # all forced boundaries is O(n_segments) CLQD calls with no benefit.
  # Detect modeNum=2 by checking if atfcut is sorted and contiguous with cutpoints.
  skip_remerge <- length(atfcut) > 0L && !is.null(atfcut) &&
    length(cutpoints) > 2L && (length(atfcut) >= length(cutpoints) - 2L)
  if (length(atfcut) && !skip_remerge) {
    newLDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2L)
    consec <- 0L; remerge_count <- 0L
    for (i in seq_len(nrow(done) - 1L)) {
      eb <- done[i,]; nb <- done[i+1L,]
      # Fast O(log n) check using sorted atfcut instead of O(n*m) intersect()
      fi <- findInterval(eb[2L], atfcut)
      has_forced_cut <- fi < length(atfcut) && atfcut[fi + 1L] <= nb[1L]
      if (has_forced_cut) {
        consec <- consec + 1L
        if (consec > 1L) { eb <- newLDblocks[remerge_count, , drop=FALSE] }
        rng <- range(c(eb, nb))
        bv  <- CLQD(geno[,rng[1L]:rng[2L]], SNPinfo[rng[1L]:rng[2L],],
                    adjN[,rng[1L]:rng[2L]],
                    CLQcut=CLQcut, clstgap=clstgap, CLQmode=CLQmode,
                    codechange=FALSE, checkLargest=checkLargest, split=split,
                    digits=digits, n_threads=n_threads, verbose=FALSE)
        bins <- seq_len(max(bv[!is.na(bv)]))
        cl   <- lapply(bins, function(x) sort(which(bv == x)))
        cl   <- Filter(length, cl)  # drop empty bins (Louvain IDs non-contiguous)
        if (length(cl)) cl <- cl[order(vapply(cl, min, numeric(1L)))]
        tmp  <- constructLDblock(cl, SNPinfo[rng[1L]:rng[2L],])
        tmp  <- tmp + (rng[1L]-1L); tmp <- tmp[order(tmp[,1L]),,drop=FALSE]
        newLDblocks[(remerge_count+1L):(remerge_count+nrow(tmp)), ] <- tmp
        remerge_count <- remerge_count + nrow(tmp)
      } else {
        consec <- 0L
        remerge_count <- remerge_count + 1L; newLDblocks[remerge_count,] <- eb
        if (i == nrow(done)-1L) { remerge_count <- remerge_count + 1L; newLDblocks[remerge_count,] <- nb }
      }
    }
    LDblocks <- newLDblocks[seq_len(remerge_count), , drop=FALSE]
  } else { LDblocks <- done }

  LDblocks <- LDblocks[!is.na(LDblocks[,1L]),,drop=FALSE]
  LDblocks <- LDblocks[order(LDblocks[,1L]),,drop=FALSE]

  # -- Resolve overlapping blocks (LD-informed split) ------------------------
  # Overlaps arise at sub-segment seams: a boundary SNP can be assigned to
  # a block in segment A and also to a block in segment B, producing two
  # blocks whose SNP index ranges overlap.
  #
  # OLD behaviour: blind union merge -- take [min(all), max(all)].
  #   Problem: merges blocks with low inter-block LD, creating an
  #   artificially large block with a deflated He and a high count of
  #   low-frequency haplotype alleles.
  #
  # NEW behaviour: LD-informed split --
  #   For each disputed SNP in the overlap zone, compute its mean r2
  #   with a sample of SNPs from the LEFT core (A-exclusive region) and
  #   from the RIGHT core (B-exclusive region). Assign it to whichever
  #   core it is in stronger LD with. The block boundary is set at the
  #   last left-assigned SNP, preserving contiguity of both blocks.
  #
  #   Fallback to union merge when either core is empty (one block fully
  #   contained in the other), where there is no LD anchor to compare.
  #
  # Parameters:
  #   k_rep = max representatives sampled from each core (bounds cost)
  #   Uses col_r2_cpp() -- O(n_ind) per call, negligible overhead.

  # -- LD-informed overlap resolution ----------------------------------------
  # Replaces the blind union merge [min(starts), max(ends)] with a per-SNP
  # LD-based assignment of disputed SNPs to one of the two overlapping blocks.
  #
  # GEOMETRY (after sorting LDblocks by start index):
  #   Block A: [sA ......... sB-1 | sB ......... eA]
  #   Block B:              [sB ......... eA | eA+1 ......... eB]
  #                          ^--- overlap zone ---^
  #   left_core:  [sA .. sB-1]   A-exclusive (always non-empty after sort)
  #   overlap:    [sB .. eA]     disputed SNPs

  # -- LD-informed overlap resolution (C++ accelerated) -------------------
  # resolve_overlap_cpp() implements the identical cumulative-score split
  # rule as the original R .resolve_overlap(), but:
  #   1. Runs entirely in C++ (no R interpreter overhead per block pair)
  #   2. Computes r2 ONLY against left_reps union right_reps (at most 2*k_rep
  #      columns), not all p columns - O(n * 2k_rep) vs O(n * p) per SNP
  #   3. Uses a lazy column cache: standardises each column at most once
  #   4. Single pass (not twice): one call is sufficient because the C++
  #      implementation processes all overlapping pairs in order; the second
  #      R call was needed to catch overlaps introduced by the first pass's
  #      boundary adjustments, which the C++ version handles within the same
  #      while loop by re-checking position i after each split.
  if (nrow(LDblocks) > 1L) {
    LDblocks <- resolve_overlap_cpp(LDblocks, adjN, k_rep = 10L, singleton_as_block = singleton_as_block)
  }

  # -- LD-informed gap closing -------------------------------------------------
  # Any inter-block gap that contains SNPs is assigned to the flanking LD blocks.
  # This prevents unassigned SNPs between adjacent blocks. Empty physical gaps
  # without SNPs remain acceptable because no SNP index exists there.
  if (isTRUE(close_gaps_with_snps) && nrow(LDblocks) > 1L) {
    LDblocks <- .close_interblock_gaps(
      blocks = LDblocks,
      adj_mat = adjN,
      k_rep = gap_k_rep
    )
  }

  # -- Handle SNPs not assigned to any clique ----------------------------------
  #
  # During CLQD clustering, some SNPs may remain NA in bv.
  # These SNPs were collected into singleton_idx_list.
  #
  # Final behaviour:
  #
  #   singleton_as_block = TRUE
  #     -> report genuine isolated SNPs as singleton LD blocks.
  #
  #   singleton_as_block = FALSE
  #     -> absorb them into neighbouring LD blocks before output construction.

  unassigned_idx <- integer(0L)

  if (exists("singleton_idx_list") &&
      length(singleton_idx_list) > 0L) {

    unassigned_idx <- unlist(
      singleton_idx_list,
      use.names = FALSE
    )

    unassigned_idx <- unique(
      as.integer(
        unassigned_idx[!is.na(unassigned_idx)]
      )
    )
  }

  if (length(unassigned_idx) > 0L && nrow(LDblocks) > 0L) {

    covered_idx <- unique(
      unlist(
        Map(
          seq.int,
          LDblocks[, 1L],
          LDblocks[, 2L]
        ),
        use.names = FALSE
      )
    )

    unassigned_idx <- setdiff(
      unassigned_idx,
      covered_idx
    )
  }

  if (length(unassigned_idx) > 0L) {

    if (isTRUE(singleton_as_block)) {

      singleton_blocks <- cbind(
        start = unassigned_idx,
        end   = unassigned_idx
      )

      LDblocks <- rbind(
        LDblocks,
        singleton_blocks
      )

    } else {

      LDblocks <- .absorb_unassigned_snps(
        blocks = LDblocks,
        unassigned_idx = unassigned_idx,
        adj_mat = adjN,
        k_rep = gap_k_rep
      )
    }

    LDblocks <- LDblocks[
      order(LDblocks[, 1L], LDblocks[, 2L]),
      ,
      drop = FALSE
    ]

    if (nrow(LDblocks) > 1L) {
      LDblocks <- resolve_overlap_cpp(
        LDblocks,
        adjN,
        k_rep = gap_k_rep,
        singleton_as_block = singleton_as_block
      )
    }

    if (!isTRUE(singleton_as_block)) {
      LDblocks <- LDblocks[
        LDblocks[, 1L] < LDblocks[, 2L],
        ,
        drop = FALSE
      ]
    }
  }

  # -- Build output data.frame -----------------------------------------------

  if (nrow(LDblocks) == 0L) {
    return(
      data.frame(
        start = integer(),
        end = integer(),
        start.rsID = character(),
        end.rsID = character(),
        start.bp = numeric(),
        end.bp = numeric(),
        n_snps = integer()
      )
    )
  }
  out <- data.frame(
    start      = LDblocks[,1L],
    end        = LDblocks[,2L],
    start.rsID = as.character(SNPinfo[[1L]][LDblocks[,1L]]),
    end.rsID   = as.character(SNPinfo[[1L]][LDblocks[,2L]]),
    start.bp   = as.numeric(SNPinfo[[2L]][LDblocks[,1L]]),
    end.bp     = as.numeric(SNPinfo[[2L]][LDblocks[,2L]])
  )
  # n_snps is computed BEFORE re-indexing while indices still map 1:1 to
  # the filtered SNPinfo rows (every index from start to end IS a SNP that
  # passed MAF filtering). After re-indexing over the full set including
  # monomorphics, end - start + 1 would over-count; computing here avoids that.
  out$n_snps <- out$end - out$start + 1L

  # Re-index over full set including monomorphic SNPs
  full_index <- stats::setNames(
    seq_len(nrow(OSNPinfo)),
    as.character(OSNPinfo[[1L]])
  )

  out$start <- unname(
    full_index[as.character(out$start.rsID)]
  )

  out$end <- unname(
    full_index[as.character(out$end.rsID)]
  )

  out <- out[
    !is.na(out$start) &
      !is.na(out$end),
    ,
    drop = FALSE
  ]

  out <- out[order(out$start),]
  if (nrow(out) > 0L) {
    # Vectorised pmin/pmax replaces t(apply(..., sort)) which is O(n) R-loop
    # and extremely slow at hundreds-of-thousands of blocks (WGS scale).
    s_tmp        <- out[["start"]];    e_tmp        <- out[["end"]]
    out[["start"]] <- pmin(s_tmp, e_tmp); out[["end"]] <- pmax(s_tmp, e_tmp)
    s_tmp        <- out[["start.bp"]]; e_tmp        <- out[["end.bp"]]
    out[["start.bp"]] <- pmin(s_tmp, e_tmp); out[["end.bp"]] <- pmax(s_tmp, e_tmp)
  }

  if (isTRUE(appendrare)) {
    prep_full <- list(V_inv_sqrt = V_inv_sqrt)
    out <- appendSGTs(out, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest,
                      CLQmode, prep_full, split, digits, n_threads)
  }

  # ---------------------------------------------------------------------------
  # Final containment cleanup
  # ---------------------------------------------------------------------------
  out <- out[
    order(out$start, out$end),
    ,
    drop = FALSE
  ]

  if (nrow(out) > 1L) {

    keep <- rep(TRUE, nrow(out))

    for (i in seq_len(nrow(out))) {

      if (!keep[i]) next

      contained <- which(
        keep &
          seq_len(nrow(out)) != i &
          out$start >= out$start[i] &
          out$end <= out$end[i]
      )

      if (length(contained)) {
        keep[contained] <- FALSE
      }
    }

    out <- out[
      keep,
      ,
      drop = FALSE
    ]
  }

  out
}
