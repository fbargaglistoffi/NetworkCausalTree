set.seed(13256)

library(NetworkCausalTree)

dataset <- data_generator_direct_indirect(
  N = 2000,
  M = 4,
  k = 80,
  p = rep(0.2, 2000),
  het = TRUE,
  h = 2,
  method_networks = "er",
  param_er = 0.1
)

result <- NetworkCausalTree(
  X = dataset[["X"]],
  Y = dataset[["Y"]],
  W = dataset[["W"]],
  A = dataset[["A"]],
  K = dataset[["K"]],
  p = dataset[["p"]],
  effect_weights = c(1, 0, 0, 0),
  ratio_disc = 0.5,
  depth = 3,
  minsize = 5,
  method = "singular",
  output = "estimation"
)

title <- expression(paste("CAUSAL TREE TARGETED TO ", tau, "(1,0;0,0)"))
cov_names <- colnames(dataset[["X"]])

plot_NCT(
  result,
  cov_names,
  title,
  output = "estimation"
)
