# Contributing to NetworkCausalTree

Thank you for your interest in contributing to NetworkCausalTree! We welcome contributions from the community and appreciate your efforts to improve this package.

## Getting Started

NetworkCausalTree is implemented in R and follows standard R development conventions. Before contributing, please:

1. Read our [Code of Conduct](CODE_OF_CONDUCT.md)
2. Check the [GitHub Issues](https://github.com/fbargaglistoffi/NetworkCausalTree/issues) page for existing issues or feature requests
3. Review the package documentation and vignettes to understand the current functionality

## Types of Contributions

We welcome various types of contributions:

### Bug Reports

If you find a bug, please report it by:

1. Opening a new issue on our [GitHub Issues](https://github.com/fbargaglistoffi/NetworkCausalTree/issues) page
2. Including a clear title and description
3. Providing a minimal reproducible example (reprex)
4. Specifying your R version, operating system, and package version
5. Describing the expected vs. actual behavior

### Feature Requests

Have an idea for a new feature? We'd love to hear it! Please:

1. Check if a similar request already exists in the Issues
2. Open a new issue with the "enhancement" label
3. Clearly describe the feature and its potential use cases
4. Explain how it would benefit users of the package

### Code Contributions

#### Setting Up Your Development Environment

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

# Clone the repository
git clone https://github.com/fbargaglistoffi/NetworkCausalTree.git
cd NetworkCausalTree

# Load the package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

#### Contribution Workflow

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/NetworkCausalTree.git
   ```
3. **Create a new branch** for your feature or bug fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```
4. **Make your changes** following our coding standards (see below)
5. **Write or update tests** to cover your changes
6. **Update documentation** as needed (including roxygen2 comments)
7. **Run checks** to ensure everything passes:
   ```r
   devtools::check()
   devtools::test()
   ```
8. **Commit your changes** with clear, descriptive commit messages:
   ```bash
   git commit -m "Add feature: brief description of what you did"
   ```
9. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```
10. **Submit a Pull Request** on GitHub with a clear description of your changes

### Documentation Improvements

Documentation improvements are always welcome! This includes:

- Fixing typos or clarifying existing documentation
- Adding examples to function documentation
- Improving or creating vignettes
- Updating the README

To update function documentation:

1. Edit the roxygen2 comments in the R source files
2. Run `devtools::document()` to regenerate the documentation
3. Preview changes with `?function_name`

### Vignettes and Examples

We encourage contributions that improve reproducibility:

- Add example scripts demonstrating specific use cases
- Create or improve vignettes showing real-world applications
- Provide datasets that illustrate network causal inference scenarios

## Coding Standards

To maintain consistency across the codebase, please follow these guidelines:

### R Style Guide

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable and function names
- Keep functions focused on a single task
- Use `<-` for assignment, not `=`
- Indent with 2 spaces (no tabs)

### Function Documentation

All exported functions must include roxygen2 documentation with:

- `@title` and `@description`
- `@param` for each parameter
- `@return` describing the return value
- `@examples` showing usage
- `@export` if the function should be user-facing

Example:
```r
#' @title Estimate Network Causal Effects
#' @description Estimates treatment effects under network interference
#' @param X Matrix of covariates
#' @param Y Outcome vector
#' @param W Treatment assignment vector
#' @return A list containing estimated effects and standard errors
#' @examples
#' result <- estimate_effects(X, Y, W)
#' @export
estimate_effects <- function(X, Y, W) {
  # Function implementation
}
```

### Testing

- All new features should include tests using `testthat`
- Tests should be placed in the `tests/testthat/` directory
- Aim for good coverage of edge cases and error conditions
- Use `expect_equal()`, `expect_error()`, etc. for assertions

Example test structure:
```r
test_that("estimate_effects works correctly", {
  # Setup test data
  X <- matrix(rnorm(100), ncol = 2)
  Y <- rnorm(50)
  W <- rbinom(50, 1, 0.5)
  
  # Run function
  result <- estimate_effects(X, Y, W)
  
  # Check results
  expect_type(result, "list")
  expect_true("effects" %in% names(result))
})
```

## Development Priorities

We are particularly interested in contributions related to:

1. **Continuous Exposures**: Extending the package to handle multiple dosage levels (e.g., vaccine dosage)
2. **Ensemble Methods**: Integration with random forest estimators for ensemble-based causal discovery (similar to grf building on causalTree)
3. **Diagnostic Tools**: Developing sensitivity analysis tools under partial interference
4. **CNI Assumption Testing**: Tools to assess whether the Conditional Network Independence assumption is reasonable for a given dataset
5. **Performance Optimization**: Improving computational efficiency for large networks
6. **Visualization Tools**: Enhanced plotting and visualization capabilities

## Questions and Support

If you have questions or need support:

1. Check the package vignettes and documentation first
2. Search existing [GitHub Issues](https://github.com/fbargaglistoffi/NetworkCausalTree/issues)
3. If your question hasn't been addressed, open a new issue with the "question" label

## Recognition

Contributors will be acknowledged in:

- The package DESCRIPTION file (for significant contributions)
- Release notes
- The project README

## License

By contributing to NetworkCausalTree, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).

## Thank You!

Your contributions help make NetworkCausalTree better for everyone. We appreciate your time and effort in improving this package!
