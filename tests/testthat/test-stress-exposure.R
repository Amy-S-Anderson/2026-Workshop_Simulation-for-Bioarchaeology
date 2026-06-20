# Unit tests for stress exposure sampling

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE)),
    lesion = rep(c(0,1), n/2)
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("stress exposure sampling raises error when annual_exposure > 1", {
  state <- create_test_state()

  expect_error(
    sample_exposure(
      pop = state,
      annual_exposure = 2
    )
  )
})


