
set.seed(13256)

library(NetworkCausalTree)

dataset <- data_generator_direct_indirect(
  N = 4000,
  M = 4,
  k = 80,
  p = rep(0.2, 4000),
  het = TRUE,
  h = 3,
  method_networks = "sf"
)

result <- NetworkCausalTree(
  X = dataset[["X"]],
  Y = dataset[["Y"]],
  W = dataset[["W"]],
  A = dataset[["A"]],
  K = dataset[["K"]],
  p = dataset[["p"]],
  effect_weights = c(0.25, 0.25, 0.25, 0.25),
  ratio_disc = 0.5,
  depth = 2,
  minsize = 5,
  method = "composite",
  output = "detection"
)

title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
cov_names <- colnames(dataset[["X"]])

plot_NCT(
  NCT = result,
  cov_names = cov_names,
  title = title,
  output = "detection"
)
