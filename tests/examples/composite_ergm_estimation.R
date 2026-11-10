set.seed(13256)

library(NetworkCausalTree)

dataset <- data_generator_direct_indirect(
  N = 1000,
  M = 4,
  k = 100,
  p = rep(0.2, 1000),
  het = TRUE,
  h = 4,
  method_networks = "ergm",
  coef_ergm = c(-0.10, 0.1),
  var_homophily_ergm = "x1"
)

result <- NetworkCausalTree(
  X = dataset[["X"]],
  Y = dataset[["Y"]],
  W = dataset[["W"]],
  A = dataset[["A"]],
  K = dataset[["K"]],
  p = dataset[["p"]],
  effect_weights = c(0.5, 0, 0, 0.5),
  ratio_disc = 0.5,
  depth = 3,
  minsize = 5,
  method = "composite",
  output = "estimation"
)

title <- expression(
  paste("CAUSAL TREE TARGETED TO ",
        tau, "(1,0;0,0) & ",
        tau, "(0,1;0,0)")
)

cov_names <- colnames(dataset[["X"]])

plot_NCT(
  NCT = result,
  cov_names = cov_names,
  title = title,
  output = "estimation"
)
