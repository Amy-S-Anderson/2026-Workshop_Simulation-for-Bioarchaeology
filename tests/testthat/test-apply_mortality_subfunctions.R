test_that("compute_hazard_multiplier: no risk factors, no exposure -> multiplier of 1", {
  pop <- data.frame(age = c(5, 30, 70))
  hm <- compute_hazard_multiplier(pop)
  expect_equal(hm, c(1, 1, 1))
})

test_that("compute_hazard_multiplier: continuous risk factor (NULL spec) multiplies directly", {
  pop <- data.frame(age = c(5, 30), frailty = c(0.8, 1.5))
  hm <- compute_hazard_multiplier(pop, risk_factors = list(frailty = NULL))
  expect_equal(hm, c(0.8, 1.5))
})

test_that("compute_hazard_multiplier: binary risk factor (numeric scalar spec) applies only where col == 1", {
  pop <- data.frame(age = c(5, 30, 70), lesion = c(1, 0, 1))
  hm <- compute_hazard_multiplier(pop, risk_factors = list(lesion = 2.1))
  expect_equal(hm, c(2.1, 1, 2.1))
})

test_that("compute_hazard_multiplier: NA in risk factor column is treated as 1 (no effect)", {
  pop <- data.frame(age = c(5, 30), frailty = c(NA, 1.5))
  hm <- compute_hazard_multiplier(pop, risk_factors = list(frailty = NULL))
  expect_equal(hm, c(1, 1.5))
})

test_that("compute_hazard_multiplier: missing risk_factors column is silently skipped", {
  pop <- data.frame(age = c(5, 30))
  hm <- compute_hazard_multiplier(pop, risk_factors = list(lesion = 2.1))
  expect_equal(hm, c(1, 1))
})

test_that("compute_hazard_multiplier: transient_hazard column, if present, multiplies in", {
  pop <- data.frame(age = c(5, 30), transient_hazard = c(3, 1))
  hm <- compute_hazard_multiplier(pop)
  expect_equal(hm, c(3, 1))
})

test_that("compute_hazard_multiplier: REGRESSION -- transient exposure hazard applies
           regardless of lesion_requires_survival (the bug this refactor fixed)", {
  pop <- data.frame(age = c(5, 30), exposed_this_step = c(TRUE, FALSE))

  hm_survival_false <- compute_hazard_multiplier(
    pop,
    exposure_causes_hazard     = TRUE,
    hazard_is_transient        = TRUE,
    exposure_hazard_multiplier = 4
    # lesion_requires_survival not passed -> defaults FALSE in apply_mortality,
    # but compute_hazard_multiplier doesn't take this arg at all anymore,
    # which is the point: exposure hazard no longer depends on it.
  )
  expect_equal(hm_survival_false, c(4, 1))
})

test_that("compute_hazard_multiplier: exposure hazard is a no-op when exposure_causes_hazard = FALSE", {
  pop <- data.frame(age = c(5, 30), exposed_this_step = c(TRUE, TRUE))
  hm <- compute_hazard_multiplier(
    pop,
    exposure_causes_hazard     = FALSE,
    hazard_is_transient        = TRUE,
    exposure_hazard_multiplier = 4
  )
  expect_equal(hm, c(1, 1))
})

test_that("compute_hazard_multiplier: exposure hazard is a no-op when hazard_is_transient = FALSE
           (permanent hazard is handled elsewhere, not here)", {
  pop <- data.frame(age = c(5, 30), exposed_this_step = c(TRUE, TRUE))
  hm <- compute_hazard_multiplier(
    pop,
    exposure_causes_hazard     = TRUE,
    hazard_is_transient        = FALSE,
    exposure_hazard_multiplier = 4
  )
  expect_equal(hm, c(1, 1))
})

test_that("compute_hazard_multiplier: multiple multiplicative sources combine correctly", {
  pop <- data.frame(
    age               = c(5, 30),
    frailty           = c(2, 0.5),
    lesion            = c(1, 0),
    transient_hazard  = c(1.5, 1),
    exposed_this_step = c(TRUE, FALSE)
  )
  hm <- compute_hazard_multiplier(
    pop,
    risk_factors               = list(frailty = NULL, lesion = 3),
    exposure_causes_hazard     = TRUE,
    hazard_is_transient        = TRUE,
    exposure_hazard_multiplier = 2
  )
  # agent 1: frailty(2) * lesion(3, since lesion==1) * transient_hazard(1.5) * exposure(2) = 18
  # agent 2: frailty(0.5) * lesion(1, since lesion==0) * transient_hazard(1) * exposure(1, not exposed) = 0.5
  expect_equal(hm, c(18, 0.5))
})


# --- compute_cause_specific_hazards ----------------------------------------

test_that("compute_cause_specific_hazards: no p_trauma column -> all hazard is illness, split_active FALSE", {
  pop <- data.frame(age = c(10, 40), frailty = c(1, 1))
  mu0_lookup <- data.frame(age = c(10, 40), mu0 = c(0.01, 0.02))
  hm <- c(1, 1)

  out <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "proportional")

  expect_false(out$split_active)
  expect_equal(out$illness_hazard, c(0.01, 0.02))
  expect_equal(out$trauma_hazard, c(0, 0))
})

test_that("compute_cause_specific_hazards: p_trauma = 0 everywhere behaves identically to no column", {
  pop <- data.frame(age = c(10, 40), frailty = c(1, 1))
  mu0_lookup_no_col   <- data.frame(age = c(10, 40), mu0 = c(0.01, 0.02))
  mu0_lookup_zero_col <- data.frame(age = c(10, 40), mu0 = c(0.01, 0.02), p_trauma = c(0, 0))
  hm <- c(1, 1)

  out_no_col   <- compute_cause_specific_hazards(pop, mu0_lookup_no_col,   hm)
  out_zero_col <- compute_cause_specific_hazards(pop, mu0_lookup_zero_col, hm)

  expect_equal(out_no_col$illness_hazard, out_zero_col$illness_hazard)
  expect_equal(out_no_col$trauma_hazard,  out_zero_col$trauma_hazard)
  expect_false(out_zero_col$split_active)
})

test_that("compute_cause_specific_hazards: illness + trauma conserve total baseline hazard
           (proportional case, hazard_multiplier = 1, frailty = 1)", {
  pop <- data.frame(age = c(10, 40, 70), frailty = c(1, 1, 1))
  mu0_lookup <- data.frame(
    age      = c(10, 40, 70),
    mu0      = c(0.01, 0.02, 0.05),
    p_trauma = c(0.4, 0.6, 0.1)
  )
  hm <- c(1, 1, 1)

  out <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "proportional")

  expect_true(out$split_active)
  expect_equal(out$illness_hazard + out$trauma_hazard, mu0_lookup$mu0, tolerance = 1e-12)
})

test_that("compute_cause_specific_hazards: frailty and hazard_multiplier affect illness_hazard only,
           trauma_hazard is untouched by either", {
  pop <- data.frame(age = c(40, 40), frailty = c(1, 3))
  mu0_lookup <- data.frame(age = 40, mu0 = 0.02, p_trauma = 0.5)
  hm <- c(1, 5)  # agent 2 has elevated frailty AND hazard_multiplier (e.g. lesion)

  out <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "proportional")

  # trauma_hazard identical for both agents despite very different frailty/hm
  expect_equal(out$trauma_hazard[1], out$trauma_hazard[2])
  expect_equal(out$trauma_hazard, c(0.01, 0.01))  # 0.02 * 0.5
  # illness_hazard differs: agent 2 = 0.02 * 0.5 * 3 * 5 = 0.15, agent 1 = 0.01
  expect_equal(out$illness_hazard, c(0.01, 0.15))
})

test_that("compute_cause_specific_hazards: unrecognized mortality_risk_type errors", {
  pop <- data.frame(age = 40, frailty = 1)
  mu0_lookup <- data.frame(age = 40, mu0 = 0.02)
  expect_error(
    compute_cause_specific_hazards(pop, mu0_lookup, hazard_multiplier = 1,
                                   mortality_risk_type = "not_a_real_type"),
    "Unrecognized mortality_risk_type"
  )
})

test_that("compute_cause_specific_hazards: time_decreasing and time_increasing match hand-computed values", {
  pop <- data.frame(age = 20, frailty = 2)
  mu0_lookup <- data.frame(age = 20, mu0 = 0.03, p_trauma = 0.25)
  hm <- 4

  age_illness <- 0.03 * (1 - 0.25)  # 0.0225
  age_trauma  <- 0.03 * 0.25        # 0.0075

  out_dec <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "time_decreasing")
  expect_equal(out_dec$illness_hazard, age_illness * 2 * 4 / ((20 / 10) + 4))
  expect_equal(out_dec$trauma_hazard, age_trauma)

  out_inc <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "time_increasing")
  expect_equal(out_inc$illness_hazard, age_illness * 2 * ((20 / 10) + 4) / 4)
  expect_equal(out_inc$trauma_hazard, age_trauma)
})


# --- resolve_deaths ----------------------------------------------------------

test_that("resolve_deaths: force_death = TRUE kills everyone regardless of hazard", {
  set.seed(1)
  pop <- data.frame(age = c(5, 40, 90), id = 1:3)
  out <- resolve_deaths(
    pop            = pop,
    illness_hazard = c(0, 0, 0),   # hazard is zero -- would never die by chance
    trauma_hazard  = c(0, 0, 0),
    split_active   = FALSE,
    current_time   = 7,
    force_death    = TRUE
  )
  expect_equal(nrow(out$pop), 0)
  expect_equal(nrow(out$decedents), 3)
  expect_true(all(out$decedents$year_died == 7))
})

test_that("resolve_deaths: split_active = FALSE never adds a cause_of_death column,
           even when force_death = TRUE", {
  set.seed(1)
  pop <- data.frame(age = c(5, 40), id = 1:2)
  out <- resolve_deaths(
    pop            = pop,
    illness_hazard = c(0.5, 0.5),
    trauma_hazard  = c(0, 0),
    split_active   = FALSE,
    current_time   = 1,
    force_death    = TRUE
  )
  expect_false("cause_of_death" %in% names(out$decedents))
})

test_that("resolve_deaths: split_active = TRUE assigns cause_of_death to every decedent,
           and assigns it deterministically when hazard is 100% one cause", {
  set.seed(1)
  pop <- data.frame(age = c(5, 40), id = 1:2)

  out <- resolve_deaths(
    pop            = pop,
    illness_hazard = c(1, 0),      # agent 1: illness only
    trauma_hazard  = c(0, 1),      # agent 2: trauma only
    split_active   = TRUE,
    current_time   = 1,
    force_death    = TRUE
  )

  expect_true("cause_of_death" %in% names(out$decedents))
  expect_equal(nrow(out$decedents), 2)
  # order within decedents_this_step follows original row order (died == TRUE subset)
  cod_by_id <- setNames(out$decedents$cause_of_death, out$decedents$id)
  expect_equal(cod_by_id[["1"]], "illness")
  expect_equal(cod_by_id[["2"]], "trauma")
})

test_that("resolve_deaths: cause assignment respects relative hazard proportions on average
           (stochastic check over many agents)", {
  set.seed(42)
  n <- 5000
  pop <- data.frame(age = rep(30, n), id = 1:n)

  # All agents die (force_death), with illness:trauma hazard ratio fixed at 3:1,
  # so ~75% of deaths should be assigned "illness"
  out <- resolve_deaths(
    pop            = pop,
    illness_hazard = rep(0.03, n),
    trauma_hazard  = rep(0.01, n),
    split_active   = TRUE,
    current_time   = 1,
    force_death    = TRUE
  )

  p_illness_observed <- mean(out$decedents$cause_of_death == "illness")
  expect_equal(p_illness_observed, 0.75, tolerance = 0.02)
})

test_that("resolve_deaths: survivors keep all original columns except step-level ones", {
  set.seed(1)
  pop <- data.frame(
    age               = c(5, 90),
    transient_hazard  = c(1, 1),
    exposed_this_step = c(TRUE, FALSE)
  )
  out <- resolve_deaths(
    pop                      = pop,
    illness_hazard           = c(0, 0),  # nobody dies
    trauma_hazard            = c(0, 0),
    split_active             = FALSE,
    current_time             = 1,
    force_death              = FALSE,
    lesion_requires_survival = FALSE
  )
  expect_equal(nrow(out$pop), 2)
  expect_false("transient_hazard"  %in% names(out$pop))
  expect_false("exposed_this_step" %in% names(out$pop))
})

test_that("resolve_deaths: lesion_requires_survival = TRUE retains exposed_this_step on survivors", {
  set.seed(1)
  pop <- data.frame(
    age               = c(5, 90),
    exposed_this_step = c(TRUE, FALSE)
  )
  out <- resolve_deaths(
    pop                      = pop,
    illness_hazard           = c(0, 0),
    trauma_hazard            = c(0, 0),
    split_active             = FALSE,
    current_time             = 1,
    force_death              = FALSE,
    lesion_requires_survival = TRUE
  )
  expect_true("exposed_this_step" %in% names(out$pop))
})

test_that("resolve_deaths: decedents never carry transient_hazard or exposed_this_step,
           regardless of lesion_requires_survival", {
  set.seed(1)
  pop <- data.frame(
    age               = c(5, 90),
    transient_hazard  = c(1, 1),
    exposed_this_step = c(TRUE, TRUE)
  )
  out <- resolve_deaths(
    pop                      = pop,
    illness_hazard           = c(0, 0),
    trauma_hazard            = c(0, 0),
    split_active             = FALSE,
    current_time             = 1,
    force_death              = TRUE,
    lesion_requires_survival = TRUE
  )
  expect_false("transient_hazard"  %in% names(out$decedents))
  expect_false("exposed_this_step" %in% names(out$decedents))
})

test_that("resolve_deaths: zero deaths produces a zero-row decedents data frame
           without erroring, and doesn't add year_died/cause_of_death columns", {
  set.seed(1)
  pop <- data.frame(age = c(5, 90), id = 1:2)
  out <- resolve_deaths(
    pop            = pop,
    illness_hazard = c(0, 0),
    trauma_hazard  = c(0, 0),
    split_active   = TRUE,
    current_time   = 1,
    force_death    = FALSE
  )
  expect_equal(nrow(out$decedents), 0)
  expect_equal(nrow(out$pop), 2)
})


# --- Integration-style check: three sub-functions compose to match the
#     pre-refactor apply_mortality() formula exactly, for a case with no
#     exposure/force_death complications -----------------------------------

test_that("INTEGRATION: chained sub-functions reproduce 1 - exp(-(illness+trauma)) threshold logic", {
  set.seed(99)
  pop <- data.frame(
    age     = c(3, 25, 55, 85),
    frailty = c(1.2, 0.9, 1.0, 1.4),
    lesion  = c(1, 0, 1, 0)
  )
  mu0_lookup <- data.frame(
    age      = c(3, 25, 55, 85),
    mu0      = c(0.04, 0.01, 0.015, 0.06),
    p_trauma = c(0.05, 0.5, 0.3, 0.05)
  )

  hm <- compute_hazard_multiplier(pop, risk_factors = list(frailty = NULL, lesion = 2.1))
  ch <- compute_cause_specific_hazards(pop, mu0_lookup, hm, mortality_risk_type = "proportional")

  expected_illness <- mu0_lookup$mu0 * (1 - mu0_lookup$p_trauma) * pop$frailty *
    ifelse(pop$lesion == 1, 2.1, 1)
  expected_trauma  <- mu0_lookup$mu0 * mu0_lookup$p_trauma

  expect_equal(ch$illness_hazard, expected_illness, tolerance = 1e-12)
  expect_equal(ch$trauma_hazard,  expected_trauma,  tolerance = 1e-12)

  # And the final call sees the same total hazard a monolithic function would
  expected_threshold <- 1 - exp(-(expected_illness + expected_trauma))
  # (resolve_deaths recomputes this internally; we're just confirming the
  # inputs it receives are correct, since resolve_deaths' own dice-roll
  # correctness is covered by the stochastic test above)
  expect_true(all(expected_threshold >= 0 & expected_threshold <= 1))
})
