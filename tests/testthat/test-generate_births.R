# Unit tests for generating new agent births

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------


# Hold population configuration parameters here for ease of reference in functions that follow. 
pop_config <- list(
  model_lesions              = FALSE,
  annual_exposure            = NULL,
  model_frailty              = FALSE,
  gammafrailty_shape         = NULL,
  # nothing below here should get called by the fertility functions. 
  gammafrailty_scale         = NULL,
  lesion_formation_window    = c(0,0),
  exposure_causes_hazard     = FALSE,
  hazard_is_transient        = FALSE,
  lesion_requires_survival   = FALSE,
  exposure_hazard_multiplier = 1
)


create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE)),
    year_born = 0
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("generate_births  raises error when n_births is NULL", {
  state <- create_test_state()
  current_time = 1
  births = NULL

  expect_error(
    generate_births(
      pop = state,
      current_time = current_time,
      pop_config = pop_config,
      n_births = births
    )
  )
})

