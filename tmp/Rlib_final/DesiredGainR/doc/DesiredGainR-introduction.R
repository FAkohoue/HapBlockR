## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(DesiredGainR)
set.seed(42)

## ----data---------------------------------------------------------------------
traits <- c("yield", "protein", "disease")
n <- 60
candidates <- data.frame(
  GenoID = paste0("G", seq_len(n)),
  Family = rep(paste0("F", 1:6), each = 10),
  yield = rnorm(n),
  protein = rnorm(n),
  disease = rnorm(n)
)
G <- cov(candidates[traits])
dimnames(G) <- list(traits, traits)

## ----dgsi---------------------------------------------------------------------
dgsi <- run_dgsi(
  init_data = candidates[c("GenoID", "Family")],
  cand_data = candidates,
  trait_cols = traits,
  dg = c(yield = 0.6, protein = 0.3, disease = 0.4),
  G = G,
  lower_is_better = "disease",
  n_select = 10,
  n_iter = 100,
  n_rep = 5,
  seed = 42
)

dgsi$best_replicate
dgsi$replicate_diagnostics
head(dgsi$ranked_geno)

## ----thresholds---------------------------------------------------------------
threshold_result <- run_dgsi(
  init_data = candidates[c("GenoID", "Family")],
  cand_data = candidates,
  trait_cols = traits,
  dg = c(yield = 0.6, protein = 0.3, disease = 0.4),
  G = G,
  lower_is_better = "disease",
  select_mode = "eligible_top_n",
  trait_min = c(yield = -0.8, protein = -0.8, disease = -0.8),
  n_select = 10,
  n_iter = 100,
  n_rep = 5,
  seed = 42
)
threshold_result$eligibility

## ----qgsi---------------------------------------------------------------------
W <- matrix(
  c(
    0.10, 0.02, -0.01,
    0.02, 0.05,  0.00,
   -0.01, 0.00,  0.08
  ),
  3,
  dimnames = list(traits, traits)
)

qgsi <- run_qgsi(
  init_data = candidates[c("GenoID", "Family")],
  gebv_data = candidates,
  trait_cols = traits,
  linear_weights = c(yield = 1, protein = 0.5, disease = 0.7),
  W = W,
  lower_is_better = "disease"
)

head(qgsi$ranked_geno)
head(qgsi$linear_contributions)
head(qgsi$quadratic_contributions)

## ----session------------------------------------------------------------------
packageVersion("DesiredGainR")
sessionInfo()

