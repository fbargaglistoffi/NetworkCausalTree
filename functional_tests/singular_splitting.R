dataset <- data_generator(N = 2000,
                          m = 40,
                          p = rep(0.2,2000),
                          het = TRUE,
                          taui = 2,
                          K = 4,
                          method_networks = "er",
                          param_er = 0.1)

result <- NetworkCausalTrees(effweights = c(1,0,0,0),
                             A = dataset[["adj_matrix"]],
                             p = dataset[["p"]],
                             fracpredictors = 1,
                             W = dataset[["W"]],
                             Y = dataset[["Y"]],
                             M = dataset[["M"]],
                             G =  dataset[["G"]],
                             X =  dataset[["X"]],
                             mdisc = 25,
                             mest = 15,
                             minpopfrac = 1,
                             depth = 3,
                             minsize = 5,
                             n_trees = 1,
                             method = "singular",
                             output = "estimation")

title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
cov_names <- colnames(dataset[["X"]])
plot(result, cov_names, title)
