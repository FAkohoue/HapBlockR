library(DesiredGainR)
library(data.table)

ext <- system.file("extdata", package = "DesiredGainR")
values <- fread(file.path(ext, "example_pheno.csv"))
traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
dg <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)
G <- stats::cov(as.matrix(values[, ..traits]))

result <- run_dgsi(
  init_data = values[, .(GenoID, Family)],
  cand_data = values,
  trait_cols = traits,
  dg = dg,
  G = G,
  lower_is_better = c("BL", "NBL", "VHB"),
  n_select = 10,
  n_iter = 200,
  n_rep = 5,
  seed = 42
)

print(result$replicate_diagnostics)
print(head(result$ranked_geno))
