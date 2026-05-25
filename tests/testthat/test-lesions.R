# Unit tests for lesion module
#
# HOW TEST FILES WORK IN R (testthat overview)
# -----------------------------------------------
# This file is discovered and run automatically by testthat when you call
# devtools::test() from your package root. Every file in tests/testthat/ that
# starts with "test_" is executed.
#
# The fundamental building block is test_that():
#
#   test_that("plain-English description of what should be true", {
#     ... setup code ...
#     expect_*(...)   # one or more expectations
#   })
#
# If every expect_*() inside a test_that() block passes, the test is green.
# If any one fails, the whole test_that() block is marked as failed, and
# testthat reports which expectation broke and what values it actually saw.
#
# Common expectations used here:
#   expect_error(expr)          — expr must throw an error
#   expect_no_error(expr)       — expr must NOT throw an error
#   expect_true(x)              — x must be TRUE
#   expect_false(x)             — x must be FALSE
#   expect_equal(x, y)          — x and y must be identical (within tolerance)
#   expect_s3_class(x, "cls")   — x must have S3 class "cls"
#   expect_type(x, "type")      — typeof(x) must equal "type"
#   expect_gt(x, y)             — x must be greater than y


# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------
# Fixtures are reusable setup objects. Defining them once at the top means
# every test_that() block can use them without repeating boilerplate.
# create_test_cohort() is a *function* rather than a plain object so that
# each test gets a fresh, independent copy of the cohort — mutations in one
# test cannot bleed into another.

#' Build a minimal cohort data frame for testing
#' @param n Number of agents
#' @param all_alive If TRUE, no agents are marked dead (default TRUE)
create_test_cohort <- function(n = 100, all_alive = TRUE) {
  data.frame(
    agent_id  = 1:n,
    age       = as.numeric(sample(0:80, n, replace = TRUE)),
    lesion    = rep(0L, n),           # integer 0/1, matching production code
    dead      = rep(!all_alive, n),
    in_sample = rep(TRUE, n)
  )
}

# Convenience wrapper: returns the Alive index vector for a cohort.
# form_lesions() takes an explicit Alive argument so tests must supply it.
alive_idx <- function(cohort) which(!cohort$dead)


# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------
# These tests check that the function refuses obviously wrong inputs before
# doing any computation. Fast failure with a clear error is better than silent
# bad output. Each expect_error() asserts that an error is thrown — the test
# PASSES when the error occurs.

test_that("lesion_formation_rate must be between 0 and 1", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 0,
                 formation_window_closes = 80,
                 lesion_formation_rate   = -0.1,
                 dx              = 1)
  )
  
  expect_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 0,
                 formation_window_closes = 80,
                 lesion_formation_rate   = 1.5,
                 dx              = 1)
  )
})

test_that("formation_window_opens must be non-negative", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = -1,
                 formation_window_closes = 10,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("formation_window_closes must be >= formation_window_opens", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 10,
                 formation_window_closes = 5,   # closes before it opens
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("Alive index must be a non-empty integer vector", {
  cohort <- create_test_cohort()
  
  expect_error(
    form_lesions(cohort, Alive = integer(0),   # empty — no living agents
                 formation_window_opens  = 0,
                 formation_window_closes = 80,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})


# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------
# These are the mirror image of the validation tests: we confirm that valid
# inputs do NOT throw errors. expect_no_error() passes when the expression
# runs cleanly to completion.

test_that("function executes with a full-lifetime window", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_no_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 0,
                 formation_window_closes = 80,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("function executes with a window that opens mid-life", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_no_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 18,
                 formation_window_closes = 80,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("function executes with a window that closes early", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_no_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 0,
                 formation_window_closes = 6,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("function executes with a bounded interval window", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  expect_no_error(
    form_lesions(cohort, Alive,
                 formation_window_opens  = 2,
                 formation_window_closes = 6,
                 lesion_formation_rate   = 0.05,
                 dx              = 1)
  )
})

test_that("function returns a data frame", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.05,
                         dx              = 1)
  
  expect_s3_class(result, "data.frame")
})


# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------
# These confirm that the function returns a data frame with the same shape and
# column types as the input. The lesion module must not add, remove, or rename
# columns, and must not change the number of rows.

test_that("output contains the same columns as input", {
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.05,
                         dx              = 1)
  
  expect_equal(sort(names(result)), sort(names(cohort)))
})

test_that("output has the same dimensions as input", {
  cohort <- create_test_cohort(n = 150)
  Alive  <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.05,
                         dx              = 1)
  
  expect_equal(dim(result), dim(cohort))
})

test_that("output columns have expected types", {
  # This test pins down the column types of the returned data frame so that a
  # refactor cannot silently change integer lesion to logical, etc.
  cohort <- create_test_cohort()
  Alive  <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.05,
                         dx              = 1)
  
  expect_type(result$agent_id,  "integer")
  expect_type(result$age,       "double")
  expect_type(result$lesion,    "integer")   # 0/1 integer, not logical
  expect_type(result$dead,      "logical")
  expect_type(result$in_sample, "logical")
})


# -----------------------------------------------------------------------------
# Lesion acquisition logic tests
# -----------------------------------------------------------------------------
# These are the most important tests: they check the *behaviour* of the
# function, not just its shape. They use controlled inputs (known ages, known
# lesion states, sometimes set.seed() for reproducibility) so the expected
# output can be stated precisely.

test_that("agents with existing lesions retain lesion status", {
  # Once a lesion is acquired it can never be lost. pmax() in form_lesions
  # guarantees this, but we verify it explicitly.
  cohort        <- create_test_cohort(n = 100)
  cohort$lesion <- rep(c(0L, 1L), 50)   # alternating: half already have lesions
  Alive         <- alive_idx(cohort)
  had_lesion    <- cohort$lesion == 1L
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)  # near-certain acquisition
  
  expect_true(all(result$lesion[had_lesion] == 1L))
})

test_that("dead agents are not processed for lesion acquisition", {
  # Agents in the Alive vector are those for whom form_lesions rolls dice.
  # Dead agents are simply absent from Alive, so their lesion column must not
  # change regardless of their age.
  cohort             <- create_test_cohort(n = 100)
  cohort$dead[1:20]  <- TRUE
  cohort$lesion[1:20]<- 0L
  cohort$age[1:20]   <- 5     # within any reasonable window
  Alive              <- alive_idx(cohort)  # excludes rows 1:20
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_true(all(result$lesion[1:20] == 0L))
})

test_that("no transitions occur when lesion_formation_rate is zero", {
  cohort        <- create_test_cohort(n = 100)
  cohort$age    <- rep(5, 100)   # everyone in a typical window
  cohort$lesion <- rep(0L, 100)
  Alive         <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0,
                         dx              = 1)
  
  expect_true(all(result$lesion == 0L))
})

test_that("near-certain rate causes most agents in window to acquire lesions", {
  # This is a stochastic test: with rate = 0.99 and n = 1000 the probability
  # of fewer than 95% converting is astronomically small, so it is safe to
  # assert without set.seed().
  cohort        <- create_test_cohort(n = 1000)
  cohort$age    <- rep(5, 1000)
  cohort$lesion <- rep(0L, 1000)
  Alive         <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_gt(mean(result$lesion), 0.95)
})


# -----------------------------------------------------------------------------
# Age window tests
# -----------------------------------------------------------------------------

test_that("agents outside the window (above closes) do not acquire lesions", {
  cohort <- data.frame(
    agent_id  = 1:100,
    age       = rep(10, 100),   # above formation_window_closes = 6
    lesion    = rep(0L, 100),
    dead      = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 6,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_true(all(result$lesion == 0L))
})

test_that("agents outside the window (below opens) do not acquire lesions", {
  cohort <- data.frame(
    agent_id  = 1:100,
    age       = rep(1, 100),    # below formation_window_opens = 5
    lesion    = rep(0L, 100),
    dead      = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 5,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_true(all(result$lesion == 0L))
})

test_that("agents inside a bounded interval window can acquire lesions", {
  cohort <- data.frame(
    agent_id  = 1:150,
    age       = c(rep(1, 50), rep(4, 50), rep(10, 50)),  # below / in / above
    lesion    = rep(0L, 150),
    dead      = rep(FALSE, 150),
    in_sample = rep(TRUE, 150)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 2,
                         formation_window_closes = 6,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  # Below window — no transitions
  expect_true(all(result$lesion[1:50]    == 0L))
  # Above window — no transitions
  expect_true(all(result$lesion[101:150] == 0L))
  # Inside window — most should transition
  expect_gt(mean(result$lesion[51:100]), 0.95)
})


# -----------------------------------------------------------------------------
# Window boundary tests
# -----------------------------------------------------------------------------
# These pin down *inclusive* vs *exclusive* behaviour at the exact boundary
# ages, matching the >= / <= conditions in form_lesions.

test_that("agents exactly at formation_window_opens can acquire lesions", {
  set.seed(42)
  cohort <- data.frame(
    agent_id  = 1:1000,
    age       = rep(18, 1000),   # exactly at opens
    lesion    = rep(0L, 1000),
    dead      = rep(FALSE, 1000),
    in_sample = rep(TRUE, 1000)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 18,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_gt(mean(result$lesion), 0.95)  # lower bound is inclusive
})

test_that("agents just below formation_window_opens cannot acquire lesions", {
  cohort <- data.frame(
    agent_id  = 1:100,
    age       = rep(17.9, 100),  # just below opens = 18
    lesion    = rep(0L, 100),
    dead      = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 18,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_true(all(result$lesion == 0L))
})

test_that("agents exactly at formation_window_closes can acquire lesions", {
  # form_lesions uses <=, so the closes boundary is *inclusive*
  set.seed(42)
  cohort <- data.frame(
    agent_id  = 1:1000,
    age       = rep(6, 1000),    # exactly at closes
    lesion    = rep(0L, 1000),
    dead      = rep(FALSE, 1000),
    in_sample = rep(TRUE, 1000)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 6,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_gt(mean(result$lesion), 0.95)  # upper bound is inclusive
})

test_that("agents just above formation_window_closes cannot acquire lesions", {
  cohort <- data.frame(
    agent_id  = 1:100,
    age       = rep(7, 100),     # just above closes = 6
    lesion    = rep(0L, 100),
    dead      = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  Alive <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 6,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_true(all(result$lesion == 0L))
})


# -----------------------------------------------------------------------------
# Age unchanged tests
# -----------------------------------------------------------------------------

test_that("agent ages are not modified by the lesion module", {
  # form_lesions must be a pure lesion-status update — it must not touch
  # any other column. We verify age specifically because an earlier version
  # of the module accidentally incremented age as a side-effect.
  cohort       <- create_test_cohort(n = 100)
  original_age <- cohort$age
  Alive        <- alive_idx(cohort)
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 80,
                         lesion_formation_rate   = 0.05,
                         dx              = 1)
  
  expect_equal(result$age, original_age)
})


# -----------------------------------------------------------------------------
# Combined eligibility tests
# -----------------------------------------------------------------------------
# This is a single, precise scenario that exercises all four combinations of
# (dead/alive) x (in-window/out-of-window) at once.

test_that("only alive agents without lesions in window can acquire lesions", {
  set.seed(42)
  cohort <- data.frame(
    agent_id  = 1:8,
    age       = as.numeric(c(3, 3, 3, 3, 10, 10, 3, 3)),
    lesion    = c(0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L),
    dead      = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
    in_sample = rep(TRUE, 8)
  )
  # Agent 1: alive, no lesion, age 3 (in [0,6])   -> CAN transition
  # Agent 2: dead,  no lesion, age 3 (in window)  -> CANNOT (dead)
  # Agent 3: alive, has lesion, age 3             -> already 1, stays 1
  # Agent 4: dead,  has lesion, age 3             -> CANNOT (dead)
  # Agent 5: alive, no lesion, age 10 (outside)  -> CANNOT (outside window)
  # Agent 6: dead,  no lesion, age 10             -> CANNOT (dead + outside)
  # Agent 7: alive, no lesion, age 3              -> CAN transition
  # Agent 8: alive, no lesion, age 3              -> CAN transition
  
  Alive  <- alive_idx(cohort)  # rows 1, 3, 5, 7, 8
  
  result <- form_lesions(cohort, Alive,
                         formation_window_opens  = 0,
                         formation_window_closes = 6,
                         lesion_formation_rate   = 0.99,
                         dx              = 1)
  
  expect_false(result$lesion[2] == 1L)   # dead agent: no change
  expect_true( result$lesion[3] == 1L)   # pre-existing lesion: retained
  expect_true( result$lesion[4] == 1L)   # dead + pre-existing: retained
  expect_false(result$lesion[5] == 1L)   # outside window: no change
  expect_false(result$lesion[6] == 1L)   # dead + outside: no change
})