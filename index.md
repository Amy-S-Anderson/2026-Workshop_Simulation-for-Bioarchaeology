# Persephone ABM

Agent-based model for simulating skeletal lesion formation and mortality
in archaeological populations, based on the Usher3 illness-death model.

## Installation

Install from GitHub using `devtools`:

``` r

# Install devtools if needed
install.packages("devtools")

# Install persephone (and demohaz dependency)
devtools::install_github("Amy-S-Anderson/persephone")
```

## Dependencies

- R (\>= 3.5)
- demohaz (installed automatically from GitHub)

## Development

### Building the Package

``` r

# Generate documentation from roxygen2 comments
devtools::document()

# Build and install locally
devtools::install()

# Check package for CRAN compliance (optional)
devtools::check()
```

### Running Tests

``` r

# Run all tests
devtools::test()

# Run a specific test file
testthat::test_file("tests/testthat/test-estimation-error.R")
```

### Loading for Development

``` r

# Load package without installing (for interactive development)
devtools::load_all()
```

## Usage

``` r

library(persephone)

# Run a simulation
result <- Simulate_Cemetery(
  cohort_size = 500,
  lesion_formation_rate = 0.10,
  formation_window_closes = 5,
  mortality_regime = CoaleDemenyWestF5
)

# Access individual outcomes
head(result$individual_outcomes)

# Access survivor data
head(result$survivors)
```
