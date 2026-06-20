# Unit tests for age increase in surviving agents at end of time step

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE))
  )
}


# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

age_pop <- function(pop) {
  pop$age <- pop$age + 1
  pop
}

test_that("age_pop raises error when age is NA", {
  state <- create_test_state()
  state$age <- NA
  
  expect_error(
    age_pop(
      pop = state
    )
  )
})



test_that("age_pop raises error when input is not a dataframe", {
  state <- create_test_state()
  state <- state$age

  expect_error(
    age_pop(
      pop = state
    )
  )
})

test_that("age_pop raises error when input dataframe does contain a column called `age`", {
  state <- create_test_state()
  state$Age <- state$age
  
  expect_error(
    age_pop(
      pop = state
    )
  )
})



test_that("age_pop raises error when age is not numeric", {
  state <- create_test_state()
  state$age <- as.character(state$age)
  
  expect_error(
    age_pop(
      pop = state
    )
  )
})


# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("age of surviving agents increases by 1", {
  state <- create_test_state()
  result <- age_pop(state)
  aged_state <- state$age + 1

  expect_equal(result$age, aged_state)
})


# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("output state contains same columns as input state", {
  state <- create_test_state()
  
  result <- age_pop(
    pop = state
  )
  expect_equal(sort(names(result)), sort(names(state)))
})



test_that("output state has same dimensions as input state", {
  state <- create_test_state(n = 150)

  result <- age_pop(
    pop = state
  )
  
  expect_equal(nrow(result), nrow(state))
  expect_equal(ncol(result), ncol(state))
  expect_equal(dim(result), dim(state))
})

