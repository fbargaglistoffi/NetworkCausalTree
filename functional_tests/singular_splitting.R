dataset <- data_generator(N = 2000,
                          K = 4,
                          m = 40,
                          p = rep(0.2,2000),
                          het = TRUE,
                          taui = 2,
                          method_networks = "er",
                          param_er = 0.1)

result <- NetworkCausalTrees(X =  dataset[["X"]],
                             Y = dataset[["Y"]],
                             W = dataset[["W"]],
                             effweights <- c(1,0,0,0),
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
                             method = "singular",
                             output = "estimation")


title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)",sep=""))
cov_names <- colnames(dataset[["X"]])
plot_NCT(result, cov_names, title)
