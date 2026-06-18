# Persephone ABM

An agent-based model (ABM) for simulating skeletal assemblage data in bioarchaeology. Individuals in a birth cohort or age-structured population face annual risks of age-specific mortality following a Siler hazard model and, optionally, annual risks of skeletal lesion formation.  
  
  Lesion status can modify mortality risk, enabling exploration of the osteological paradox and evaluation of survival analysis methods or multistate illness-death models applied to cemetery data.


## Installation

Install from GitHub using `devtools`:

```r
# Install devtools if needed
install.packages("devtools")

# Install persephone (and demohaz dependency)
devtools::install_github("Amy-S-Anderson/persephone")
```

## Dependencies

- R (>= 3.5)
- demohaz (installed automatically from GitHub)

## Development

### Building the Package

```r
# Generate documentation from roxygen2 comments
devtools::document()

# Build and install locally
devtools::install()

# Check package for CRAN compliance (optional)
devtools::check()
```

### Running Tests

```r
# Run all tests
devtools::test()

# Run a specific test file
testthat::test_file("tests/testthat/test-estimation-error.R")
```

### Loading for Development

```r
# Load package without installing (for interactive development)
devtools::load_all()
```

## Usage
The main function in the package is `Simulate_Cemetery()`. 
```r
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
head(result$annual_census)
```
