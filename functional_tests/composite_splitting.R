set.seed(1)

dataset <- data_generator(N = 2000,
                          K = 4,
                          m = 40,
                          p = rep(0.2,2000),
                          het = FALSE,
                          h = 0,
                          method_networks = "sf")

result <- NetworkCausalTrees(X =  dataset[["X"]],
                             Y = dataset[["Y"]],
                             W = dataset[["W"]],
                             A = dataset[["A"]],
                             M = dataset[["M"]],
                             p = dataset[["p"]],
                             effweights = c(0.25,0.25,0.25,0.25),
                             ratio_dis = 0.5,
                             depth = 3,
                             minsize = 5,
                             method = "composite",
                             output = "estimation")

title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
cov_names <- colnames(dataset[["X"]])
plot_NCT(result, cov_names, title)
