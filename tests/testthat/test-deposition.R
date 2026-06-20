# Unit tests for deposition module

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

# Deposition parameters: cutoff age
deposition_param_cutoff <- c(2.5)  # cutoff at age 2.5

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

test_that("deposition module raises error when deposition_model is NULL", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_deposition(
      decedents = state,
      deposition_model = NULL,
      deposition_param = deposition_param_cutoff,
      dx = dx
    )
  )
})

test_that("deposition module raises error when deposition_model is unsupported", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_deposition(
      decedents = state,
      deposition_model = "unknown",
      deposition_param = deposition_param_cutoff,
      dx = dx
    )
  )

  expect_error(
    apply_deposition(
      decedents = state,
      deposition_model = "siler",
      deposition_param = deposition_param_cutoff,
      dx = dx
    )
  )
})

test_that("cutoff model requires length-1 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_deposition(
      decedents = state,
      deposition_model = "cutoff",
      deposition_param = c(1, 2),  # wrong length
      dx = dx
    )
  )

  expect_error(
    apply_deposition(
      decedents = state,
      deposition_model = "cutoff",
      deposition_param = c(),  # empty
      dx = dx
    )
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("cutoff model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_deposition(
      decedents = state,
      deposition_model = "cutoff",
      deposition_param = deposition_param_cutoff,
      dx = dx
    )
  )
})

test_that("deposition module returns a dataframe", {
  state <- create_test_state()
  dx <- 1

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param_cutoff,
    dx = dx
  )

  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("output state has expected column types", {
  state <- create_test_state()
  dx <- 1

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param_cutoff,
    dx = dx
  )

  expect_type(result$agent_id, "integer")
  expect_type(result$age, "double")
  expect_type(result$was_deposited, "logical")
})

# -----------------------------------------------------------------------------
# Deposition logic tests
# -----------------------------------------------------------------------------

test_that("dead agents below cutoff are not deposited", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(c(0, 1, 2, 2.4, 2.5, 3, 4, 5, 6, 7)),
    lesion = rep(FALSE, 10)
  )
  dx <- 1
  deposition_param <- c(2.5)  # cutoff at age 2.5

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param,
    dx = dx
  )

  # Agents at ages 0, 1, 2, 2.4 (below cutoff) should not be deposited
  expect_false(result$was_deposited[1])
  expect_false(result$was_deposited[2])
  expect_false(result$was_deposited[3])
  expect_false(result$was_deposited[4])
})



test_that("agent exactly at cutoff age is deposited", {
  state <- data.frame(
    agent_id = 1:5,
    age = as.numeric(c(2.49, 2.5, 2.51, 4, 5))
  )
  dx <- 1
  deposition_param <- c(2.5)  # cutoff at age 2.5

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param,
    dx = dx
  )

  # Age 2.49 is below cutoff
  expect_false(result$was_deposited[1])
  # Age 2.5 is exactly at cutoff - should be deposited
  expect_true(result$was_deposited[2])
  # Ages 2.51, 4, 5 are above cutoff
  expect_true(result$was_deposited[3])
  expect_true(result$was_deposited[4])
  expect_true(result$was_deposited[5])
})

test_that("cutoff of 1 also works correctly", {
  state <- data.frame(
    agent_id = 1:5,
    age = as.numeric(c(0.5, 1.0, 1.5, 2, 3)),
    lesion = rep(FALSE, 5)
  )
  dx <- 1
  deposition_param <- c(1)  # cutoff at age 1

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param,
    dx = dx
  )

  # Age 0.5 is below cutoff
  expect_false(result$was_deposited[1])
  # Ages 1.0, 1.5, 2, 3 are at or above cutoff
  expect_true(result$was_deposited[2])
  expect_true(result$was_deposited[3])
  expect_true(result$was_deposited[4])
  expect_true(result$was_deposited[5])
})

test_that("cutoff logic is deterministic (no randomness)", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(c(2, 3.5), 50)),
    lesion = rep(FALSE, 100)
  )
  dx <- 1
  deposition_param <- c(2.5)

  # Run twice with different seeds
  set.seed(111)
  result1 <- apply_deposition(state, "cutoff", deposition_param, dx)
  set.seed(999)
  result2 <- apply_deposition(state, "cutoff", deposition_param, dx)

  # Results should be identical (deterministic)
  expect_equal(result1$was_deposited, result2$was_deposited)
})

# -----------------------------------------------------------------------------
# State preservation tests
# -----------------------------------------------------------------------------

test_that("agent ages are not modified by deposition module", {
  state <- create_test_state(n = 100)
  dx <- 1

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param_cutoff,
    dx = dx
  )

  expect_equal(result$age, state$age)
})

test_that("lesion status is not modified by deposition module", {
  state <- create_test_state(n = 100)
  dx <- 1

  result <- apply_deposition(
    decedents = state,
    deposition_model = "cutoff",
    deposition_param = deposition_param_cutoff,
    dx = dx
  )

  expect_equal(result$lesion, state$lesion)
})

