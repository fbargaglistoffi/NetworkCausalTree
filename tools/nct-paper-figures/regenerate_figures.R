# Regenerate Figures 2 and 3 from the final JOSS manuscript examples.
# The data-generating and fitting arguments match paper 14.pdf, pages 6-8.
# Added reproducibility seeds are fixed in advance, not selected by outcomes.
# Run at repository root: Rscript tools/nct-paper-figures/regenerate_figures.R
# Optional argument: direct or composite (otherwise both examples are run).
# No package R functions are modified by this script.

options(warn = 1)
source_commit <- "c93e2c3fa1e8b0d3d28d4fdbefb154c036ed749d"
if (dir.exists("R") && file.exists("DESCRIPTION") &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(NetworkCausalTree)
}

args <- commandArgs(trailingOnly = TRUE)
examples <- if (length(args)) args else c("direct", "composite")
stopifnot(all(examples %in% c("direct", "composite")))
output_root <- Sys.getenv("NCT_FIGURE_OUTPUT", "nct-paper-figures")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

run_example <- function(example) {
  direct <- identical(example, "direct")
  figure <- if (direct) 2L else 3L
  data_seed <- if (direct) 20261116L else 20261117L
  split_seed <- if (direct) 20261118L else 20261119L
  out <- file.path(output_root, paste0("figure", figure))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  warnings_seen <- character()
  started <- Sys.time()
  cat("BEGIN EXAMPLE:", example, "\n")
  cat("Source commit:", source_commit, "\n")
  cat("Data seed:", data_seed, "; split seed:", split_seed, "\n")

  set.seed(data_seed)
  generator <- if (direct) data_generator_direct else data_generator_direct_indirect
  dataset <- generator(
    N = 4000, M = 4, k = 80, p = rep(0.2, 4000),
    het = TRUE, h = 2, method_networks = "er", param_er = 0.1
  )
  cat("Generated observations after the generator's default isolate removal:",
      length(dataset$Y), "\n")
  cat("Clusters:", length(unique(dataset$K)), "\n")
  stopifnot(length(unique(dataset$K)) == 80L,
            all(rowSums(dataset$A) > 0),
            length(dataset$Y) <= 4000L)
  saveRDS(dataset, file.path(out, "dataset.rds"), compress = "xz")

  set.seed(split_seed)
  result <- withCallingHandlers(
    NetworkCausalTree(
      X = dataset[["X"]], Y = dataset[["Y"]], W = dataset[["W"]],
      A = dataset[["A"]], K = dataset[["K"]], p = dataset[["p"]],
      effect_weights = if (direct) c(1, 0, 0, 0) else rep(0.25, 4),
      ratio_disc = 0.5, depth = 2, minsize = 5,
      method = if (direct) "singular" else "composite",
      output = "estimation"
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
    }
  )
  saveRDS(result, file.path(out, "fitted_tree.rds"))
  write.csv(result, file.path(out, "fitted_tree.csv"), row.names = FALSE)
  cat("FITTED TREE BEGIN\n")
  print(result)
  cat("FITTED TREE END\n")

  root <- which(result$NODE == 1L)
  stopifnot(length(root) == 1L)
  total <- result$NOBS_TR + result$NOBS_EST
  stopifnot(total[root] == length(dataset$Y),
            result$NOBS_TR[root] > 0, result$NOBS_EST[root] > 0)
  leaves <- result$TERMINAL == "LEAF"
  stopifnot(sum(result$NOBS_TR[leaves]) == result$NOBS_TR[root],
            sum(result$NOBS_EST[leaves]) == result$NOBS_EST[root])
  percent_fun <- get("nct_node_percentages", envir = asNamespace("NetworkCausalTree"))
  percentages <- percent_fun(result)
  stopifnot(isTRUE(all.equal(percentages, 100 * total / total[root])),
            abs(percentages[root] - 100) < 1e-12)
  effect_columns <- paste0("EFF", c("1000", "1101", "1110", "0100"), "_EST")
  se_columns <- paste0("SE", c("1000", "1101", "1110", "0100"), "_EST")
  stopifnot(all(is.finite(as.matrix(result[, effect_columns, drop = FALSE]))),
            all(is.finite(as.matrix(result[, se_columns, drop = FALSE]))))
  labels <- data.frame(
    NODE = result$NODE, FILTER = result$FILTER,
    direct_effect = result$EFF1000_EST, direct_se = result$SE1000_EST,
    spillover_effect = result$EFF0100_EST, spillover_se = result$SE0100_EST,
    N_discovery = result$NOBS_TR, N_estimation = result$NOBS_EST,
    N_total = total, percent_of_root = percentages
  )
  write.csv(labels, file.path(out, "plotted_values.csv"), row.names = FALSE)

  plot_title <- if (direct) {
    expression(paste("NETWORK CAUSAL TREE: ", tau, "(1,0;0,0)"))
  } else {
    "NETWORK CAUSAL TREE: COMPOSITE CRITERION"
  }
  draw_plot <- function() {
    plot_NCT(
      NCT = result, cov_names = colnames(dataset[["X"]]), title = plot_title,
      output = "estimation", effect_color_nodes = if (direct) "1000" else "0100",
      vertex_size = 32, vertex_size2 = 24, vertex_label_cex = 1,
      edge_label_cex = 0.95, edge_arrow_size = 0.45,
      main_cex = 1.1, digits = 4
    )
  }
  basename <- paste0("Figure_", figure, if (direct) "_direct" else "_composite")
  grDevices::cairo_pdf(file.path(out, paste0(basename, ".pdf")),
                      width = 9, height = 7, pointsize = 11, family = "sans")
  draw_plot()
  grDevices::dev.off()
  grDevices::png(file.path(out, paste0(basename, ".png")),
                 width = 9, height = 7, units = "in", res = 600,
                 pointsize = 11, type = "cairo", bg = "white")
  draw_plot()
  grDevices::dev.off()
  grDevices::svg(file.path(out, paste0(basename, ".svg")),
                 width = 9, height = 7, pointsize = 11, family = "sans")
  draw_plot()
  grDevices::dev.off()

  manifest <- list(
    source_commit = source_commit, example = example, figure = figure,
    data_seed = data_seed, split_seed = split_seed,
    requested_N = 4000, actual_N = length(dataset$Y),
    removed_isolates = 4000 - length(dataset$Y), clusters = 80,
    root_discovery_N = result$NOBS_TR[root], root_estimation_N = result$NOBS_EST[root],
    root_percentage = percentages[root], leaves = sum(leaves),
    warnings = unique(warnings_seen),
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    note = paste("Same generator and model arguments as the manuscript.",
                 "Node N combines discovery and estimation samples, as in plot_NCT.",
                 "Effect estimates and standard errors use only held-out estimation units.",
                 "The generator's default remove_isolates=TRUE was retained.")
  )
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE,
                       auto_unbox = TRUE, digits = 16)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(unique(warnings_seen), file.path(out, "warnings.txt"))
  cat("MANIFEST BEGIN\n")
  cat(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
  cat("MANIFEST END\n")
  cat("END EXAMPLE:", example, "\n")
  invisible(result)
}

invisible(lapply(examples, run_example))
