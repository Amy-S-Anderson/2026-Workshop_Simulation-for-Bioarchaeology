# ------------------------------------------------------------------------------
# Stable age distribution initialization for Persephone ABM
# ------------------------------------------------------------------------------


# --- Internal helpers (not exported) ---

#' Create the initial pop data frame
#' 
#' 
#' Compute discrete Siler survivorship l(x)
#'
#' Returns the probability of surviving from birth to each age x, under the
#' Siler hazard model. l(0) is defined as 1 (everyone is alive at birth).
#'
#' @param ages Integer vector of ages (typically 0:max_age)
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric vector of survivorship values, same length as ages
#' @keywords internal
#' @export
compute_siler_survivorship <- function(ages, mortality_regime) {
  hazards        <- compute_siler_risk(ages, mortality_regime)
  survival_probs <- 1 - hazards
  cumprod(c(1, survival_probs[-length(survival_probs)]))
}



#' Compute the stable age distribution from Siler parameters
#'
#' Derives the proportion of the living population at each age under the
#' stable age distribution, which is determined by the mortality regime and
#' the population growth rate r. For a stationary population (r = 0) this
#' reduces to the survivorship function l(x) normalized to sum to 1.
#'
#' @param ages Integer vector of ages to include (e.g. 0:80)
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @param r Numeric. Annual population growth rate. Default 0 (stationary).
#' @return Numeric vector of proportions summing to 1, same length as ages
#' @keywords internal
#' @export
compute_stable_age_distribution <- function(ages, mortality_regime, r = 0) {
  lx <- compute_siler_survivorship(ages, mortality_regime)
  
  # Weight each age class by the growth-rate discount factor e^(-rx).
  # When r = 0 all weights are 1 and this reduces to lx alone.
  weights <- exp(-r * ages) * lx
  
  # Normalize to a proper probability distribution
  weights / sum(weights)
}




#' Initialize a cohort with a stable age distribution
#'
#' Rather than starting all agents at age 0, this seeds the population with
#' agents distributed across ages according to the stable age distribution
#' implied by the mortality regime and growth rate. This allows the simulation
#' to begin in (approximate) demographic equilibrium without a burn-in period.
#'
#' Each agent's age is drawn by sampling from the stable age distribution.
#' Lesion status at initialization is set to FALSE for all agents; if you want
#' agents to enter with lesions reflecting prior exposure, a separate
#' prevalence-at-age initialization step would be needed.
#'
#' @param pop0_size Integer. Number of agents in the starting population.
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#' @param r Numeric. Annual population growth rate. Default 0 (stationary).
#' @param max_age Integer. Oldest age class to include in the distribution.
#'   Default 100.
#' @return A data frame with the same structure as create_cohort(), but with
#'   ages drawn from the stable age distribution rather than all set to 0.
#'
#' @examples
#' pop <- create_pop_stable_age(
#'   pop0_size    = 500,
#'   mortality_regime = CoaleDemenyWestF5,
#'   r              = 0.01
#' )
#'
#' @export
create_pop_stable_age <- function(pop0_size,
                                  mortality_regime,
                                  r       = 0,
                                  lesion_formation_rate = NULL,
                                  max_age = 100) {
  ages <- 0:max_age
  
  sad <- compute_stable_age_distribution(ages, mortality_regime, r)
  
  # Sample each agent's starting age from the stable age distribution
  starting_ages <- sample(ages, size = pop0_size, replace = TRUE, prob = sad)
  
  pop0 <- data.frame(
    agent_id     = 1:pop0_size,
    age          = as.numeric(starting_ages)
  )
  pop0$year_born <- if_else(pop0$age == 0, 0, NA_real_)
  
  if(!is.null(lesion_formation_rate )){
    pop0 <- pop0 %>%
      mutate(lesion = 0L) %>%
      relocate(lesion, .after = age)
  }
  return(pop0)
}

### This function is then called inside the 'create_pop() function