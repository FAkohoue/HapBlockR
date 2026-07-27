library(DGQGSI)
library(data.table)

ext <- system.file("extdata", package = "DGQGSI")
pheno <- fread(file.path(ext, "example_pheno.csv"))
gebv  <- fread(file.path(ext, "example_gebv.csv"))
traits <- c("YLD","MY","MI","BL","NBL","VHB")
dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

res <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = pheno[, .(GenoID, Family)],
  cand_data = pheno[, c("GenoID", traits), with = FALSE],
  ref_data = pheno[, c("GenoID", traits), with = FALSE],
  gebv_data = gebv,
  trait_cols = traits,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  trait_min_sd = c(YLD=0.2, MY=0.1, MI=0.1, BL=0.1, NBL=0.1, VHB=0.1),
  dg_scale_traits = FALSE,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits = FALSE,
  n_iter = 50,
  n_rep = 3,
  merge_outputs = TRUE,
  debug = FALSE
)

print(head(res$comparison_result$comparison_table))
