dataset <- data_generator(N = 2000,
                          K = 4,
                          m = 40,
                          p = rep(0.2,2000),
                          het = FALSE,
                          taui = 0,
                          method_networks = "sf")

result <- NetworkCausalTrees(X =  dataset[["X"]],
                             Y = dataset[["Y"]],
                             W = dataset[["W"]],
                             effweights <- c(0.25,0.25,0.25,0.25),
                             A = dataset[["A"]],
                             G =  dataset[["G"]],
                             M = dataset[["M"]],
                             p = dataset[["p"]],
                             mdisc = 25,
                             mest = 15,
                             minpopfrac = 1,
                             fracpredictors = 1,
                             n_trees = 1,
                             depth = 3,
                             minsize = 5,
                             method = "composite",
                             output = "estimation")

title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
cov_names <- colnames(dataset[["X"]])

plot(result, cov_names, title)
