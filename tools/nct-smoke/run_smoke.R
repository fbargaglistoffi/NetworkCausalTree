# Small diagnostic simulation, not an AOAS reproduction or power study.
# Uses the actual revised R functions loaded from the CI workspace.
options(warn = 1, width = 160)
pkgload::load_all('.', quiet = TRUE)
outdir <- 'nct-smoke-results'
dir.create(outdir, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
B <- as.integer(Sys.getenv('NCT_SMOKE_REPS', '12'))
nclusters <- 80L
cluster_size <- 4L
N <- nclusters * cluster_size
minsize <- 3L
methods <- c('singular', 'composite', 'penalized')
effects <- c('1000', '1101', '1110', '0100')
writeLines(c(paste('replications =', B), paste('N =', N),
             paste('clusters =', nclusters), 'cluster topology = four-cycle',
             'assignment probabilities = 0.25, 0.35, 0.45, 0.55 within every cluster',
             'discovery proportion = 0.5; depth = 1; minsize = 3',
             'direct effect = 1 + 3 * X.age; spillover = 0.5 + X.age',
             'Y = direct * W + spillover * G + 0.2 * X.noise + Normal(0, 0.5^2)',
             'data seed = 20261100 + replicate; split seed = 20261200 + replicate'),
           file.path(outdir, 'config.txt'))

A <- matrix(0, N, N)
for (k in seq_len(nclusters)) {
  ids <- (4L * k - 3L):(4L * k)
  for (j in seq_len(4L)) {
    a <- ids[j]; b <- ids[if (j == 4L) 1L else j + 1L]
    A[a, b] <- A[b, a] <- 1
  }
}
K <- rep(paste0('cluster_', 1000L + 7L * seq_len(nclusters)), each = 4L)
x <- rep(rep(0:1, each = nclusters / 2L), each = 4L)
p <- rep(c(0.25, 0.35, 0.45, 0.55), nclusters)
neighbors <- lapply(seq_len(N), function(i) which(A[i, ] == 1))
q0 <- vapply(neighbors, function(ids) prod(1 - p[ids]), numeric(1))
tau <- 1 + 3 * x
spill <- 0.5 + x
filter_fun <- get('evaluate_nct_filter', envir = asNamespace('NetworkCausalTree'))
fit_rows <- list(); effect_rows <- list(); group_rows <- list(); warning_rows <- list()
fit_objects <- list()

for (r in seq_len(B)) {
  set.seed(20261100L + r)
  X <- data.frame(X.age = x, X.noise = rbinom(N, 1L, 0.5))
  W <- rbinom(N, 1L, p)
  G <- as.numeric(A %*% W > 0)
  Y <- tau * W + spill * G + 0.2 * X$X.noise + rnorm(N, sd = 0.5)
  split_seed <- 20261200L + r
  set.seed(split_seed)
  disc_clusters <- sample(unique(K), round(nclusters * 0.5), replace = FALSE)
  disc <- which(K %in% disc_clusters)
  est <- which(!K %in% disc_clusters)
  canonical <- data.frame(idunit = seq_len(N), W = W, G = G, Y = Y,
                          X.1 = X$X.age, X.2 = X$X.noise)
  for (method in methods) {
    weights <- if (method == 'composite') c(0.5, 0, 0, 0.5) else c(1, 0, 0, 0)
    warns <- character()
    set.seed(split_seed)
    start <- proc.time()[['elapsed']]
    fit <- tryCatch(withCallingHandlers(
      NetworkCausalTree(X = X, Y = Y, W = W, A = A, K = K, p = p,
                       effect_weights = weights, ratio_disc = 0.5,
                       depth = 1, minsize = minsize, method = method,
                       output = 'estimation'),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart('muffleWarning') }),
      error = function(e) e)
    seconds <- proc.time()[['elapsed']] - start
    if (inherits(fit, 'error')) {
      fit_rows[[length(fit_rows) + 1L]] <- data.frame(replication = r, method = method,
        success = FALSE, error = conditionMessage(fit), leaves = NA, recovered = FALSE,
        partition_ok = FALSE, minsize_ok = FALSE, honest_counts_ok = FALSE,
        finite_estimates_and_se = FALSE, max_manual_HT_error = NA_real_,
        warnings = length(warns), seconds = seconds)
      cat(sprintf('FIT %02d %-9s ERROR %s\n', r, method, conditionMessage(fit)))
      next
    }
    leaf_ids <- which(fit$TERMINAL == 'LEAF')
    membership <- vapply(leaf_ids, function(j) {
      filter_fun(as.character(fit$FILTER[j]), canonical)
    }, logical(N))
    if (is.null(dim(membership))) membership <- matrix(membership, ncol = 1L)
    partition_ok <- all(rowSums(membership) == 1L)
    size_ok <- all(vapply(seq_along(leaf_ids), function(j) {
      ids <- disc[membership[disc, j]]
      all(table(factor(W[ids], levels = 0:1), factor(G[ids], levels = 0:1)) >= minsize)
    }, logical(1)))
    # A root-only tree is legitimate; child feasibility applies only after a split.
    if (length(leaf_ids) == 1L && fit$NODE[leaf_ids] == 1L) size_ok <- TRUE
    root <- which(fit$NODE == 1L)
    honest_ok <- identical(as.integer(fit$NOBS_TR[root]), as.integer(length(disc))) &&
      identical(as.integer(fit$NOBS_EST[root]), as.integer(length(est)))
    recovered <- length(leaf_ids) == 2L &&
      all(vapply(seq_along(leaf_ids), function(j) length(unique(x[membership[, j]])) == 1L, logical(1)))
    estimate_cols <- paste0('EFF', effects, '_EST')
    se_cols <- paste0('SE', effects, '_EST')
    finite_ok <- all(is.finite(as.matrix(fit[, c(estimate_cols, se_cols)]))) &&
      all(as.matrix(fit[, se_cols]) >= 0)
    manual_error <- 0
    predicted <- matrix(NA_real_, N, length(effects), dimnames = list(NULL, effects))
    for (j in seq_len(nrow(fit))) {
      inside <- filter_fun(as.character(fit$FILTER[j]), canonical)
      ids <- est[inside[est]]
      mu <- function(w, g) {
        prob <- (if (w == 1) p[ids] else 1 - p[ids]) *
          (if (g == 1) 1 - q0[ids] else q0[ids])
        sum((W[ids] == w & G[ids] == g) * Y[ids] / prob) / length(ids)
      }
      expected <- c(mu(1, 0) - mu(0, 0), mu(1, 1) - mu(0, 1),
                    mu(1, 1) - mu(1, 0), mu(0, 1) - mu(0, 0))
      actual <- as.numeric(fit[j, estimate_cols])
      manual_error <- max(manual_error, abs(actual - expected))
      if (j %in% leaf_ids) {
        predicted[inside, ] <- matrix(rep(actual, each = sum(inside)), ncol = length(effects))
        truth <- c(mean(tau[ids]), mean(tau[ids]), mean(spill[ids]), mean(spill[ids]))
        for (h in seq_along(effects)) {
          se <- fit[[se_cols[h]]][j]
          effect_rows[[length(effect_rows) + 1L]] <- data.frame(
            replication = r, method = method, node = fit$NODE[j], effect = effects[h],
            n_est = length(ids), truth = truth[h], estimate = actual[h], se = se,
            covered = fit[[paste0('CI', effects[h], '_LOWER')]][j] <= truth[h] &&
                      fit[[paste0('CI', effects[h], '_UPPER')]][j] >= truth[h])
        }
      }
    }
    for (g in 0:1) for (h in seq_along(effects)) {
      ids <- est[x[est] == g]
      group_rows[[length(group_rows) + 1L]] <- data.frame(replication = r, method = method,
        subgroup = g, effect = effects[h],
        truth = if (h <= 2) 1 + 3 * g else 0.5 + g,
        estimate = mean(predicted[ids, h]))
    }
    fit_rows[[length(fit_rows) + 1L]] <- data.frame(replication = r, method = method,
      success = TRUE, error = '', leaves = length(leaf_ids), recovered = recovered,
      partition_ok = partition_ok, minsize_ok = size_ok, honest_counts_ok = honest_ok,
      finite_estimates_and_se = finite_ok, max_manual_HT_error = manual_error,
      warnings = length(warns), seconds = seconds)
    if (length(warns)) warning_rows[[length(warning_rows) + 1L]] <-
      data.frame(replication = r, method = method, warning = warns)
    fit_objects[[paste(r, method, sep = '_')]] <- fit
    cat(sprintf('FIT %02d %-9s leaves=%d recovered=%s HT_error=%.3g finite=%s seconds=%.2f warnings=%d\n',
                r, method, length(leaf_ids), recovered, manual_error, finite_ok, seconds, length(warns)))
  }
}
fits <- do.call(rbind, fit_rows)
leaf_effects <- do.call(rbind, effect_rows)
groups <- do.call(rbind, group_rows)
write.csv(fits, file.path(outdir, 'fits.csv'), row.names = FALSE)
write.csv(leaf_effects, file.path(outdir, 'leaf_effects.csv'), row.names = FALSE)
write.csv(groups, file.path(outdir, 'group_estimates.csv'), row.names = FALSE)
saveRDS(fit_objects, file.path(outdir, 'fitted_trees.rds'))
if (length(warning_rows)) write.csv(do.call(rbind, warning_rows), file.path(outdir, 'warnings.csv'), row.names = FALSE)
summary_rows <- lapply(methods, function(m) {
  z <- fits[fits$method == m, ]
  data.frame(method = m, runs = nrow(z), errors = sum(!z$success),
             recovered = sum(z$recovered), invariant_failures = sum(!z$partition_ok | !z$minsize_ok |
               !z$honest_counts_ok | !z$finite_estimates_and_se),
             max_manual_HT_error = max(z$max_manual_HT_error, na.rm = TRUE),
             warning_count = sum(z$warnings), elapsed_seconds = sum(z$seconds))
})
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(outdir, 'method_summary.csv'), row.names = FALSE)
cat('\nMETHOD_SUMMARY_BEGIN\n'); print(summary, row.names = FALSE); cat('METHOD_SUMMARY_END\n')
if (!is.null(groups)) {
  group_summary <- do.call(rbind, lapply(split(groups, list(groups$method, groups$subgroup, groups$effect), drop = TRUE), function(z) {
    data.frame(method = z$method[1], subgroup = z$subgroup[1], effect = z$effect[1],
               truth = z$truth[1], mean_estimate = mean(z$estimate),
               bias = mean(z$estimate - z$truth), rmse = sqrt(mean((z$estimate - z$truth)^2)),
               mcse_mean = sd(z$estimate) / sqrt(nrow(z)), n = nrow(z))
  }))
  write.csv(group_summary, file.path(outdir, 'group_summary.csv'), row.names = FALSE)
  cat('\nGROUP_SUMMARY_BEGIN\n'); print(group_summary, row.names = FALSE); cat('GROUP_SUMMARY_END\n')
  coverage_summary <- aggregate(covered ~ method + effect, leaf_effects, mean)
  write.csv(coverage_summary, file.path(outdir, 'exploratory_coverage.csv'), row.names = FALSE)
  cat('\nEXPLORATORY_COVERAGE_BEGIN\n'); print(coverage_summary, row.names = FALSE); cat('EXPLORATORY_COVERAGE_END\n')
}
cat('\nSmall diagnostic simulation only: no comprehensive performance or coverage claim.\n')
if (any(!fits$success) || any(summary$invariant_failures) || any(summary$max_manual_HT_error > 1e-9)) {
  stop('Smoke simulation identified failures; inspect the saved results.')
}
cat('SMOKE_VALIDATION_COMPLETED\n')
