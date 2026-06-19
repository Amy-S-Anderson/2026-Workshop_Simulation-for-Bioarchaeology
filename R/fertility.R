# ------------------------------------------------------------------------------
# Fertility module for Persephone ABM
# ------------------------------------------------------------------------------


#' Compute a trapezoid age-specific fertility schedule from a TFR
#'
#' The schedule rises linearly from 0 at age_start to a peak at age_peak_start,
#' remains flat to age_peak_end, then falls linearly back to 0 at age_end.
#' The peak rate is derived analytically so that the sum of ASFRs across all
#' ages equals the TFR (i.e. the schedule integrates correctly by construction).
#'
#' @param ages Integer vector of ages to compute the schedule over
#' @param tfr Numeric. Total Fertility Rate (expected lifetime births per woman)
#' @param age_start Numeric. Age at which fertility begins. Default 15.
#' @param age_peak_start Numeric. Age at which fertility reaches its peak. Default 25.
#' @param age_peak_end Numeric. Age at which fertility begins to decline. Default 35.
#' @param age_end Numeric. Age at which fertility reaches zero. Default 45.
#' @return Named numeric vector of annual age-specific fertility rates
#' @keywords internal
#' @export
compute_trapezoid_asfr <- function(ages,
                                   tfr,
                                   age_start      = 15,
                                   age_peak_start = 25,
                                   age_peak_end   = 35,
                                   age_end        = 45) {
  # Derive peak ASFR so that the trapezoid area equals the TFR.
  # Area = peak * (flat_width + ramp_up_width/2 + ramp_down_width/2)
  flat_width      <- age_peak_end   - age_peak_start
  ramp_up_width   <- age_peak_start - age_start
  ramp_down_width <- age_end        - age_peak_end
  area_per_peak   <- flat_width + ramp_up_width / 2 + ramp_down_width / 2
  peak_asfr       <- tfr / area_per_peak
  
  asfr <- numeric(length(ages))
  
  for (i in seq_along(ages)) {
    a <- ages[i]
    asfr[i] <- dplyr::case_when(
      a <  age_start      ~ 0,
      a <  age_peak_start ~ peak_asfr * (a - age_start)      / ramp_up_width,
      a <= age_peak_end   ~ peak_asfr,
      a <  age_end        ~ peak_asfr * (age_end - a)         / ramp_down_width,
      TRUE                ~ 0
    )
  }
  
  setNames(asfr, ages)
}



#' Generate new agent rows to append to the pop
#'
#' Creates n_births new agents with age 0 and no lesion, with agent_ids
#' continuing sequentially from the current maximum in the population. All other
#' columns are initialized to match the structure of create_pop().
#'
#' @param pop Population data frame (used to determine next agent_id)
#' @param current_time Numeric. Time counter for ABM main loop. 
#' @param n_births Integer. Number of new agents to create.
#' @param pop_config A list of parameters for the initial population so that new agents will have matching trait columns. 
#' @return Data frame with n_births rows ready to rbind() onto the population
#' @keywords internal
#' @export

generate_births <- function(pop, current_time, n_births, pop_config) {
  if (n_births == 0L) return(pop[0, ])
  
  max_id <- max(pop$agent_id)
  
  new_agents <- data.frame(
    agent_id = seq(max_id + 1, max_id + n_births),
    age      = 0,
    year_born = current_time
  )
  
  if (pop_config$model_lesions) {
    new_agents$lesion <- 0
  }
  if (pop_config$model_frailty) {
    new_agents$frailty          <- rgamma(n_births,
                                          shape = pop_config$gammafrailty_shape,
                                          scale = pop_config$gammafrailty_scale)
    new_agents$acquired_frailty <- NA_real_
  }
  
  if ("n_stress_events" %in% names(pop)) {
    new_agents$n_stress_events <- 0L
  }
  
  new_agents[, names(pop), drop = FALSE]
}



#' Apply fertility to the living population for one timestep
#'
#' Determines how many births occur this timestep and appends new agents to
#' the cohort. The number of births is drawn from a Poisson distribution whose
#' rate equals the sum of age-specific fertility rates across all reproductive-
#' age women in the living population.
#'
#' Women are approximated as half of all living agents aged between age_start
#' and age_end (inclusive). This assumption can be replaced once sex is added
#' as a column to the cohort data frame, by filtering on cohort$sex == "F"
#' instead.
#'
#' Birth counts are drawn as Poisson rather than Bernoulli because multiple
#' births can occur at a given age in a population of any size, and the Poisson
#' approximation to the sum of many independent Bernoulli trials is exact in
#' the limit of large population and small per-individual rate.
#'
#' @param cohort Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param tfr Numeric. Total Fertility Rate.
#' @param asfr Named numeric vector of age-specific fertility rates, as
#'   returned by compute_trapezoid_asfr(). The names must be character
#'   representations of integer ages.
#' @param dx Numeric. Timestep size. Scales fertility rates proportionally.
#'   Default 1 (annual timestep).
#' @return Updated cohort data frame with new agents appended
#' @keywords internal
#' @export
apply_fertility <- function(cohort, Alive, tfr, asfr, dx = 1) {
  
  # Ages of all living agents
  ages_alive <- cohort$age[Alive]
  
  # Identify reproductive-age agents and approximate females as half of them
  repro_ages    <- as.numeric(names(asfr))
  in_repro_window <- ages_alive %in% repro_ages
  
  repro_agents  <- Alive[in_repro_window]
  repro_ages_actual <- cohort$age[repro_agents]
  
  # Approximate females as half the reproductive-age living population.
  # When sex becomes a cohort column, replace this block with:
  #   female_agents <- repro_agents[cohort$sex[repro_agents] == "F"]
  #   female_ages   <- cohort$age[female_agents]
  n_repro       <- length(repro_agents)
  n_female_approx <- round(n_repro / 2)
  
  if (n_female_approx == 0L) return(cohort)  # no reproductive-age women
  
  # Expected births = sum of ASFR at each female's age, scaled by dx.
  # We sample ages for our approximate females proportional to their
  # representation in the reproductive-age pool, then look up their ASFR.
  female_age_sample <- sample(repro_ages_actual,
                              size    = n_female_approx,
                              replace = FALSE)
  
  asfr_values    <- asfr[as.character(female_age_sample)]
  expected_births <- sum(asfr_values * dx)
  
  # Draw realized births from a Poisson distribution
  n_births <- rpois(1, lambda = expected_births)
  
  # Append new agents to the cohort and return
  new_agents <- generate_births(cohort, n_births)
  cohort <- rbind(cohort, new_agents)
  
  cohort
}