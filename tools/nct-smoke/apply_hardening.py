"""Apply the runtime edits from hardening change sets 01--04 in CI only.
The original bundle's documentation edits are not needed for this smoke run.
No source edits are committed or pushed by this script.
"""
from pathlib import Path
import hashlib

root = Path(__file__).resolve().parents[2]
updated = {}

def replace(path, old, new):
    target = root / path
    text = updated.get(target, target.read_text())
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one exact anchor in {path}; found {count}")
    updated[target] = text.replace(old, new)

replace('R/compute_OF_Split.R',
'''                             N_tr = length(W), N_est = NULL) {''',
'''                             N_tr = length(W), N_est = NULL,
                             minsize = 1) {

  if (!is.numeric(minsize) || length(minsize) != 1L ||
      !is.finite(minsize) || minsize < 1 || minsize != floor(minsize)) {
    stop("minsize must be a positive integer.")
  }''')
replace('R/compute_OF_Split.R',
'''      if (any(left_table < 1) || any(right_table < 1)) {''',
'''      # Enforce feasibility before selecting the best-scoring candidate.
      if (any(left_table < minsize) || any(right_table < minsize)) {''')
replace('R/identify_partitions_nct.R',
'''        N_tr = N_tr,
        N_est = N_est
      )''',
'''        N_tr = N_tr,
        N_est = N_est,
        minsize = minsize
      )''')
replace('R/NetworkCausalTree.R',
'''  N <- length(W)

  if (is.null(K)''',
'''  N <- length(W)

  if (!is.numeric(minsize) || length(minsize) != 1L ||
      !is.finite(minsize) || minsize < 1 || minsize != floor(minsize)) {
    stop("minsize must be a positive integer.")
  }

  if (is.null(K)''')

replace('R/compute_OF_Value.R',
'''    if (any(!is.finite(population_effects))) {
      stop("population_effects contains non-finite values: cannot normalize composite OF")
    }
    if (any(population_effects == 0, na.rm = TRUE)) {
      stop("population_effects contains zero: cannot normalize composite OF")
    }

    inof <- alpha * ((EffTau1000(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[1]) ^ 2) +
            beta * ((EffTau1101(N = N, W = W, G = G, Y = Y,
                                p = p, Ne_list = Ne_list,
                                unit_ids = unit_ids) ^ 2) /
                      (population_effects[2]) ^ 2) +
            gamma * ((EffTau1110(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[3]) ^ 2) +
            delta * ((EffTau0100(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[4]) ^ 2)''',
'''    weights <- c(alpha, beta, gamma, delta)
    if (any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
      stop("Composite weights must be finite, non-negative, and include a positive weight.")
    }
    active <- which(weights > 0)
    if (any(!is.finite(population_effects[active]))) {
      stop("Selected population_effects contains non-finite values: cannot normalize composite OF")
    }
    if (any(population_effects[active] == 0)) {
      stop("Selected population_effects contains zero: cannot normalize composite OF")
    }

    effect_functions <- list(EffTau1000, EffTau1101, EffTau1110, EffTau0100)
    inof <- sum(vapply(active, function(j) {
      effect <- effect_functions[[j]](
        N = N, W = W, G = G, Y = Y, p = p,
        Ne_list = Ne_list, unit_ids = unit_ids
      )
      weights[j] * (effect / population_effects[j]) ^ 2
    }, numeric(1)))''')

replace('R/sprout_nct.R',
r'''    if (!all(grepl("^X\\.", colnames(X)))) {
      colnames(X) <- paste0("X.", seq_len(ncol(X)))
    }''',
'''    # Stored filters use positional names accepted by evaluate_nct_filter().
    colnames(X) <- paste0("X.", seq_len(ncol(X)))''')
for path in ['R/identify_partitions_nct.R', 'R/compute_effects_nct.R']:
    replace(path,
r'''    cn <- colnames(X)
    if (is.null(cn) || length(cn) != ncol(X) || !all(grepl("^X\\.", cn))) {
      colnames(X) <- paste0("X.", seq_len(ncol(X)))
    }''',
'''    # Match discovery and estimation columns by their original position.
    colnames(X) <- paste0("X.", seq_len(ncol(X)))''')

replace('R/NetworkCausalTree.R',
'''  if (anyNA(A)) {
    stop("A cannot contain missing values.")
  }
''',
'''  if (anyNA(A)) {
    stop("A cannot contain missing values.")
  }

  if ((!is.numeric(A) && !is.logical(A)) ||
      any(!is.finite(A)) || any(A != 0 & A != 1)) {
    stop("A must contain only finite binary entries (0 or 1).")
  }
''')
replace('R/NetworkCausalTree.R',
'''  # Sample the clusters to be assigned to the discovery set
''',
'''  # All four contrasts are reported, so each of their exposure states must
  # be supported for every focal unit. Do not silently redefine the population.
  validate_nct_exposure_support(seq_len(N), p, Ne_list)

  # Sample the clusters to be assigned to the discovery set
''')

helper = root / 'R/validate_nct_exposure_support.R'
if helper.exists():
    raise RuntimeError('Exposure-support helper already exists')
updated[helper] = '''# Internal input-domain check for the four contrasts reported by NCT.
# A zero pairwise probability is allowed; a zero marginal probability is not.
validate_nct_exposure_support <- function(unit_ids, p, Ne_list) {
  if (any(lengths(Ne_list[unit_ids]) == 0L)) {
    stop(
      "Exposure positivity is not satisfied: isolated units cannot have G = 1. ",
      "The current NetworkCausalTree interface reports all four contrasts ",
      "for the supplied target population and therefore requires non-isolated units. ",
      "No units have been removed automatically.",
      call. = FALSE
    )
  }
  for (w in 0:1) {
    for (g in 0:1) {
      probability <- pi(unit_ids, w, g, p, Ne_list)
      if (any(!is.finite(probability)) || any(probability <= 0)) {
        stop(
          "Exposure probabilities must be finite and strictly positive for ",
          "every unit and each reported exposure state. A zero may also ",
          "indicate numerical underflow; no denominator has been regularized.",
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}
'''
for target, content in updated.items():
    target.write_text(content)
    print(hashlib.sha256(content.encode()).hexdigest(), target.relative_to(root))
print('Applied all runtime edits from hardening change sets 01--04.')
