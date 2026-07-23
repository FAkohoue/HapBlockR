## tests/testthat/test-selection-clusters.R
## -----------------------------------------------------------------------------
## Tests for cluster_selection_groups() / plot_selection_clusters()
## (R/parent_selection.R): PC-variance-threshold clustering (hierarchical or
## k-means) of individuals, cross-tabulated against an arbitrary number of
## named selection strategies (GA, TS, OCS, core-collection, ...) -- the
## higher-dimensional, more rigorous companion to
## plot_parent_selection_pca()'s PC1/PC2-only visual cross-check.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Fixtures ------------------------------------------------------------------
# Three well-separated synthetic groups (large mean offsets relative to
# within-group noise) so hierarchical/k-means clustering recovers them
# reliably -- this lets tests check that each TRUE group concentrates in a
# single dominant computed cluster, without asserting exact (arbitrary)
# cluster label identity.

set.seed(42)
.n_per_grp <- 6L
.n_ind     <- 3L * .n_per_grp
.true_grp  <- rep(c("A", "B", "C"), each = .n_per_grp)
.ids       <- paste0("ind", seq_len(.n_ind))

# feature_matrix: individuals x 8 "blocks", group means offset by 20 (large
# relative to unit noise) on the first 2 columns only, small noise elsewhere.
.centres <- matrix(0, nrow = 3L, ncol = 8L)
.centres[1, 1:2] <- c(20, 0)
.centres[2, 1:2] <- c(-20, 20)
.centres[3, 1:2] <- c(-20, -20)
.fmat <- .centres[rep(1:3, each = .n_per_grp), ] +
  matrix(rnorm(.n_ind * 8), .n_ind, 8)
dimnames(.fmat) <- list(.ids, paste0("block", 1:8))

# G: a relationship matrix built from the same clustered data (so it carries
# the same 3-group structure), forced symmetric PSD-like via tcrossprod.
.G <- tcrossprod(scale(.fmat, scale = FALSE)) / 8
diag(.G) <- diag(.G) + 1e-6
dimnames(.G) <- list(.ids, .ids)

.ga_sel <- .ids[c(1, 2, 7, 8)]           # 2 from group A, 2 from group B
.ts_sel <- .ids[c(2, 3, 13, 14)]         # overlaps ind2 with GA, 2 from group C
.ocs_sel  <- .ids[c(1, 8, 9, 15)]        # a 3rd strategy, its own overlaps
.core_sel <- .ids[c(1, 7, 13)]           # a 4th strategy, one per true group
.two_groups  <- list(GA = .ga_sel, TS = .ts_sel)
.four_groups <- list(GA = .ga_sel, TS = .ts_sel, OCS = .ocs_sel, Core = .core_sel)

# ==============================================================================
# cluster_selection_groups(): input validation
# ==============================================================================

test_that("cluster_selection_groups(): errors when neither G nor feature_matrix is supplied", {
  expect_error(
    cluster_selection_groups(groups = .two_groups, n_clusters = 3L),
    "Supply either"
  )
})

test_that("cluster_selection_groups(): errors when both G and feature_matrix are supplied", {
  expect_error(
    cluster_selection_groups(G = .G, feature_matrix = .fmat,
                             groups = .two_groups, n_clusters = 3L),
    "not both"
  )
})

test_that("cluster_selection_groups(): errors when groups is missing", {
  expect_error(
    cluster_selection_groups(G = .G, n_clusters = 3L),
    "groups"
  )
})

test_that("cluster_selection_groups(): errors when groups is an unnamed list", {
  expect_error(
    cluster_selection_groups(G = .G, groups = list(.ga_sel, .ts_sel), n_clusters = 3L),
    "NAMED list"
  )
})

test_that("cluster_selection_groups(): errors when groups has duplicate names", {
  expect_error(
    cluster_selection_groups(G = .G, groups = list(GA = .ga_sel, GA = .ts_sel),
                             n_clusters = 3L),
    "unique"
  )
})

test_that("cluster_selection_groups(): errors when a groups element is not a character vector", {
  expect_error(
    cluster_selection_groups(G = .G, groups = list(GA = .ga_sel, TS = 1:3),
                             n_clusters = 3L),
    "character vector"
  )
})

test_that("cluster_selection_groups(): errors when n_clusters is missing", {
  expect_error(
    cluster_selection_groups(G = .G, groups = .two_groups),
    "n_clusters"
  )
})

test_that("cluster_selection_groups(): errors when n_clusters < 2", {
  expect_error(
    cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 1L),
    "n_clusters"
  )
})

test_that("cluster_selection_groups(): errors when n_clusters exceeds n individuals", {
  expect_error(
    cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = .n_ind + 1L),
    "cannot exceed"
  )
})

test_that("cluster_selection_groups(): errors when variance_threshold is out of (0, 1]", {
  expect_error(
    cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 3L,
                             variance_threshold = 0),
    "variance_threshold"
  )
  expect_error(
    cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 3L,
                             variance_threshold = 1.5),
    "variance_threshold"
  )
})

test_that("cluster_selection_groups(): errors on unnamed G", {
  expect_error(
    cluster_selection_groups(G = matrix(rnorm(9), 3, 3),
                             groups = list(GA = character(0)), n_clusters = 2L),
    "row/column names"
  )
})

test_that("cluster_selection_groups(): errors on feature_matrix without row names", {
  fmat_noname <- .fmat
  rownames(fmat_noname) <- NULL
  expect_error(
    cluster_selection_groups(feature_matrix = fmat_noname, groups = .two_groups,
                             n_clusters = 3L),
    "row names"
  )
})

test_that("cluster_selection_groups(): messages about unknown IDs in a group", {
  expect_message(
    cluster_selection_groups(G = .G, groups = list(GA = c(.ga_sel, "not_a_real_id")),
                             n_clusters = 3L, verbose = TRUE),
    "not among"
  )
})

# ==============================================================================
# cluster_selection_groups(): hierarchical clustering, structural correctness
# ==============================================================================

test_that("cluster_selection_groups() (hierarchical, G, 2 groups): returns correctly-shaped output", {
  res <- cluster_selection_groups(G = .G, groups = .two_groups,
                                  variance_threshold = 0.95, method = "hierarchical",
                                  n_clusters = 3L, verbose = FALSE)
  expect_equal(res$method, "hierarchical")
  expect_equal(res$space, "population diversity space")
  expect_equal(res$groups, c("GA", "TS"))
  expect_true(is.factor(res$cluster))
  expect_equal(length(res$cluster), .n_ind)
  expect_equal(nlevels(res$cluster), 3L)
  expect_true(res$n_pcs_retained >= 2L)
  expect_true(res$variance_explained >= 0.95 - 1e-8)
  expect_s3_class(res$cluster_fit, "hclust")

  # table structural checks
  expect_equal(nrow(res$table), 3L)
  expect_equal(sum(res$table$n_total), .n_ind)
  expect_equal(sum(res$table$pct_population), 100, tolerance = 1e-8)
  expect_equal(sum(res$table$n_GA), length(.ga_sel))
  expect_equal(sum(res$table$n_TS), length(.ts_sel))
  expect_equal(res$table$prop_GA,
              ifelse(res$table$n_total > 0, res$table$n_GA / res$table$n_total, NA_real_))
})

test_that("cluster_selection_groups() (hierarchical, 4 groups): every group gets its own n_*/prop_* columns", {
  res <- cluster_selection_groups(G = .G, groups = .four_groups,
                                  method = "hierarchical", n_clusters = 3L,
                                  verbose = FALSE)
  expect_equal(res$groups, c("GA", "TS", "OCS", "Core"))
  for (nm in c("GA", "TS", "OCS", "Core")) {
    expect_true(paste0("n_", nm) %in% names(res$table))
    expect_true(paste0("prop_", nm) %in% names(res$table))
  }
  expect_equal(sum(res$table$n_OCS), length(.ocs_sel))
  expect_equal(sum(res$table$n_Core), length(.core_sel))
  # ind1 is selected by BOTH GA and OCS (.ga_sel / .ocs_sel both include
  # .ids[1]) -- it must be counted in n_GA's cluster-row AND n_OCS's
  # cluster-row independently, not forced into a single mutually exclusive
  # category the way the old 2-group Neither/Both scheme would have.
  ind1_row <- res$table[res$table$Cluster == as.character(res$cluster["ind1"]), ]
  expect_true(ind1_row$n_GA[1] >= 1L)
  expect_true(ind1_row$n_OCS[1] >= 1L)
})

test_that("cluster_selection_groups() (hierarchical): recovers the 3 well-separated true groups", {
  res <- cluster_selection_groups(G = .G, groups = .two_groups,
                                  variance_threshold = 0.95, method = "hierarchical",
                                  n_clusters = 3L, verbose = FALSE)
  # Each true group's individuals should be (almost) entirely in one
  # computed cluster -- check via cross-tab purity, not exact label identity
  # (cluster "Cluster 1" vs "Cluster 2" naming is arbitrary).
  cross <- table(true = .true_grp, computed = res$cluster)
  purity <- sum(apply(cross, 1, max)) / .n_ind
  expect_gte(purity, 0.9)
})

test_that("cluster_selection_groups(): seed is ignored (with a message) for hierarchical", {
  expect_message(
    cluster_selection_groups(G = .G, groups = .two_groups,
                             method = "hierarchical", n_clusters = 3L,
                             seed = 1L, verbose = TRUE),
    "ignored"
  )
})

# ==============================================================================
# cluster_selection_groups(): k-means clustering
# ==============================================================================

test_that("cluster_selection_groups() (kmeans, feature_matrix): returns correctly-shaped output", {
  res <- cluster_selection_groups(feature_matrix = .fmat, groups = .four_groups,
                                  variance_threshold = 0.95, method = "kmeans",
                                  n_clusters = 3L, seed = 1L, verbose = FALSE)
  expect_equal(res$method, "kmeans")
  expect_equal(res$space, "target-block feature space")
  expect_equal(nlevels(res$cluster), 3L)
  expect_s3_class(res$cluster_fit, "kmeans")
})

test_that("cluster_selection_groups() (kmeans): same seed gives reproducible clustering", {
  res1 <- cluster_selection_groups(feature_matrix = .fmat, groups = .two_groups,
                                   method = "kmeans", n_clusters = 3L, seed = 7L,
                                   verbose = FALSE)
  res2 <- cluster_selection_groups(feature_matrix = .fmat, groups = .two_groups,
                                   method = "kmeans", n_clusters = 3L, seed = 7L,
                                   verbose = FALSE)
  expect_equal(as.character(res1$cluster), as.character(res2$cluster))
})

# ==============================================================================
# cluster_selection_groups(): variance_threshold controls n_pcs_retained
# ==============================================================================

test_that("cluster_selection_groups(): a low variance_threshold hits the 2-PC floor", {
  res <- cluster_selection_groups(G = .G, groups = .two_groups,
                                  variance_threshold = 0.01, method = "hierarchical",
                                  n_clusters = 3L, verbose = FALSE)
  expect_equal(res$n_pcs_retained, 2L)
})

test_that("cluster_selection_groups(): a higher variance_threshold retains more PCs", {
  res_low  <- cluster_selection_groups(G = .G, groups = .two_groups,
                                       variance_threshold = 0.5, method = "hierarchical",
                                       n_clusters = 3L, verbose = FALSE)
  res_high <- cluster_selection_groups(G = .G, groups = .two_groups,
                                       variance_threshold = 0.99, method = "hierarchical",
                                       n_clusters = 3L, verbose = FALSE)
  expect_true(res_high$n_pcs_retained >= res_low$n_pcs_retained)
})

# ==============================================================================
# plot_selection_clusters()
# ==============================================================================

test_that("plot_selection_clusters(): runs and returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  res <- cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 3L,
                                  verbose = FALSE)
  p <- plot_selection_clusters(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_selection_clusters(): errors on malformed cluster_res", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_selection_clusters(list(foo = 1)), "cluster_selection_groups")
})

test_that("plot_selection_clusters(): errors when save_path does not end in .pdf", {
  skip_if_not_installed("ggplot2")
  res <- cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 3L,
                                  verbose = FALSE)
  expect_error(plot_selection_clusters(res, save_path = tempfile(fileext = ".png")),
              "\\.pdf")
})

test_that("plot_selection_clusters(): save_path writes a real PDF file", {
  skip_if_not_installed("ggplot2")
  res <- cluster_selection_groups(G = .G, groups = .two_groups, n_clusters = 3L,
                                  verbose = FALSE)
  out_pdf <- tempfile(fileext = ".pdf")
  plot_selection_clusters(res, save_path = out_pdf)
  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)
})
