# --------------------------------------------------------------------------------------------------------------
# PERSEPHONE: Agent-Based Model for Bioarchaeology
# --------------------------------------------------------------------------------------------------------------
#
# This file contains the core simulation engine and its helper functions.
# Individual agents:
#   - Are born as a starting pop of a specified size
#   - Face a stable annual risk of forming a skeletal lesion within a specified age range
#   - May experience different age-specific mortality depending on lesion status
#   - Age-specific mortality risks follow a Siler function


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
compute_siler_survivorship <- function(ages, mortality_regime) {
  hazards <- sapply(ages, compute_siler_risk, mortality_regime = mortality_regime)
  # Probability of surviving each discrete year, then cumulative product.
  # l(0) = 1 by definition, so we prepend 1 and drop the last value to align
  # with the input age vector.
  survival_probs <- 1 - hazards
  lx <- cumprod(c(1, survival_probs[-length(survival_probs)]))
  lx
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
    age          = as.numeric(starting_ages),
    dead         = FALSE,
    was_deposited = FALSE,
    in_sample    = TRUE
  )
  if(!is.null(lesion_formation_rate )){
    pop0 <- pop0 %>%
      mutate(lesion = 0L) %>%
      relocate(lesion, .after = age)
  }
  return(pop0)
}

#' @param pop_size Number of agents in the starting pop
#' @param age_structured Logical: is this an age-structured population? (if not, it is a single age-cohort)
#' @param lesion_formation_rate Numeric, but functions as Logical. It will call the param value for this argument from the Simulate_Cemetery function call. If lesion_formation_rate has any numeric value, then a column for lesion presence is initialized; no one has lesions at time = 0. 
#' @param r Numeric, the population growth rate
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return A data frame with columns: agent_id, age, lesion, dead, in_sample
#' @keywords internal
create_pop <- function(pop0_size, age_structured, 
                       lesion_formation_rate = NULL, 
                       r = 0, mortality_regime) {
  if(age_structured == TRUE){
    pop0 <- create_pop_stable_age(pop0_size,
                                  mortality_regime,
                                  lesion_formation_rate = lesion_formation_rate,
                                  r       = r, # pop. growth rate
                                  max_age = 100)
  }
  
  if(age_structured == FALSE){
    pop0 <- data.frame(agent_id = 1:pop0_size,
                       age = 0,
                       dead = FALSE,
                       was_deposited = FALSE,
                       in_sample = TRUE
    )
    if(!is.null(lesion_formation_rate)){
      pop0 <- pop0 %>%
        mutate(lesion = 0L) %>%
        relocate(lesion, .after = age) # change position of lesion column so it sits to the right of 'age'
    }
  }
  return(pop0)
}
#' Age up all living agents to the current timestep
#' @param pop Population data frame
#' @param k Current timestep (age to assign)
#' @return Updated pop data frame
#' @keywords internal
age_pop <- function(pop, k) {
  Alive <- which(!pop$dead)
  pop$age[Alive] <- k
  pop
}

#' Compute Siler hazard for a given age
#' @param ages Integer vector, if population is age-structured
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric hazard value
#' @keywords internal
compute_siler_risk <- function(ages, mortality_regime) {
  if(length(unique(ages)) == 1){ # if all agents are the same age
    mortality_regime$a1 * exp(-mortality_regime$b1 * ages[1]) + # calculate a scalar
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages[1])
  } else{ # otherwise, calculate a vector of age-dependent values, length = n living agents
    mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }
}


#' Vectorized lesion formation across all living agents
#'
#' Replaces the per-agent form_lesion() loop with a single vectorized draw.
#' Each living agent gets its own uniform draw; conditions are evaluated
#' element-wise so individual-level heterogeneity is preserved.
#'
#' @param pop Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param formation_window_opens Age at which lesions can start forming
#' @param formation_window_closes Age at which lesions stop forming
#' @param lesion_formation_rate Annual probability of lesion formation
#' @return Updated pop data frame
#' @keywords internal
form_lesions <- function(pop, Alive,
                             formation_window_opens,
                             formation_window_closes,
                             lesion_formation_rate) {
  # One draw per living agent
  stress <- runif(length(Alive), 0, 1)
  
  ages   <- pop$age[Alive]
  in_window <- ages >= formation_window_opens & ages <= formation_window_closes
  
  # An agent acquires a lesion if: in the formation window, stress roll passes,
  # AND it does not already have a lesion (pmax preserves previously acquired lesions)
  new_lesion <- as.integer(in_window & stress <= lesion_formation_rate)
  pop$lesion[Alive] <- pmax(pop$lesion[Alive], new_lesion)
  
  pop
}

#' Vectorized mortality across all living agents
#'
#'
#' @param pop Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param age_based_risk Baseline Siler mortality risk at current age (scalar)
#' @param mortality_risk_type One of "proportional", "time_decreasing", "time_increasing"
#' @param relative_mortality_risk Multiplier for lesion-bearing individuals
#' @return Updated pop data frame
#' @keywords internal
apply_mortality <- function(pop, Alive,
                                age_based_risk,
                                mortality_risk_type,
                                relative_mortality_risk) {
  death_dice <- runif(length(Alive), 0, 1)
  
  ages   <- pop$age[Alive]
  lesion <- pop$lesion[Alive]
  
  # Compute the effective death threshold for each agent.
  # case_when() evaluates conditions row-wise (element-wise on vectors), so
  # each agent follows exactly one branch — logically equivalent to the
  # original nested ifelse() but readable and easy to extend.
  if(length(lesion) > 0){
  threshold <- dplyr::case_when(
    lesion == 0 ~
      age_based_risk,
    
    lesion == 1 & mortality_risk_type == "proportional" ~
      age_based_risk * relative_mortality_risk,
    
    lesion == 1 & mortality_risk_type == "time_decreasing" ~
      age_based_risk * relative_mortality_risk / ((ages / 10) + relative_mortality_risk),
    
    lesion == 1 & mortality_risk_type == "time_increasing" ~
      age_based_risk * ((ages / 10) + relative_mortality_risk) / relative_mortality_risk,
    
    .default = age_based_risk   # fallback: treat as no-lesion baseline
  )
  } else{
  threshold = age_based_risk
    
  }
  
  pop$dead[Alive] <- pop$dead[Alive] | (death_dice < threshold)
  
  pop
}

#' Record survivor snapshot for the current timestep. This function is only built to handle a cohort, not a full age-structured population. *This could/should be changed*
#' @param pop Population data frame
#' @param k Current timestep
#' @param lesion_formation_rate Annual rate of lesion formation (NULL if not specified)
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal
record_cohort_survivors <- function(pop, k, lesion_formation_rate) {
  n_alive  <- sum(!pop$dead & pop$age == k)
  n_lesion <- sum(!pop$dead & pop$lesion == 1 & pop$age == k)
  if(is.null(lesion_formation_rate)){
    data.frame(Age = k,
               Alive = n_alive)
  }
  else{
    data.frame(Age = k,
               Alive = n_alive,
               Lesion = n_lesion,
               Lesion_perc = ifelse(n_alive == 0, NA, round(n_lesion / n_alive * 100, 1)))
  }
}

#' Finalize the pop into a cemetery
#'
#' Ages remaining survivors one final step and marks all as dead.
#' @param pop Population data frame
#' @param k Last completed timestep
#' @return Data frame with all agents marked dead
#' @keywords internal
finalize_cemetery <- function(pop, k) {
  k <- k + 1
  Alive <- which(!pop$dead)
  pop$age[Alive] <- k
  pop$dead <- TRUE
  pop
}


# --- Main simulation function (exported) ---

#' Simulate a cemetery using the Persephone ABM
#'
#' Runs the Persephone agent-based model: a birth pop ages through time,
#' facing annual risks of skeletal lesion formation and Siler-model mortality.
#' Individuals with lesions may experience modified mortality risk. The
#' simulation ends when fewer than 10 individuals remain alive.
#'
#' @param pop0_size Integer. Number of individuals in the starting pop.
#' @param lesion_formation_rate Numeric. Annual probability of developing a
#'   lesion (between 0 and 1).
#' @param formation_window_opens Numeric. Age at which lesions can start
#'   forming. Default 0.
#' @param formation_window_closes Numeric. Age at which new lesions stop forming.
#' @param mortality_risk_type Character. How lesions modify mortality:
#'   "proportional", "time_decreasing", or "time_increasing".
#' @param relative_mortality_risk Numeric. Mortality multiplier for individuals
#'   with lesions. 1 = no effect, 2 = double risk. Default 1.
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#' @param pop_growth_rate Numeric. The population growth rate
#' @param age_structured Logical. Is the starting population age-structured, or an age cohort?
#' @return A list with two elements:
#'   \describe{
#'     \item{individual_outcomes}{Data frame of all individuals with age at
#'       death and lesion status.}
#'     \item{survivors}{Data frame of survivor counts and lesion prevalence
#'       at each age.}
#'   }
#'
#' @examples
#' result <- Simulate_Cemetery(
#'   pop0_size = 500,
#'   lesion_formation_rate = 0.10,
#'   formation_window_closes = 5,
#'   mortality_regime = CoaleDemenyWestF5
#' )
#'
#' @export
Simulate_Cemetery <- function(pop0_size,
                              lesion_formation_rate = NULL,
                              formation_window_opens = 0,
                              formation_window_closes,
                              mortality_risk_type = "proportional",
                              relative_mortality_risk = 1,
                              mortality_regime,
                              pop_growth_rate = 0, # default to stationary population
                              age_structured = TRUE,
                              deposition_param = 0,
                              taphonomy_regime,
                              loss_strength = 'no_decay',
                              age_noise = FALSE) {
  pop <- create_pop(pop0_size, age_structured = age_structured, 
                    r = pop_growth_rate, lesion_formation_rate = lesion_formation_rate,
                    mortality_regime = mortality_regime)
  
  k <- 0  # Initialize time counter
  
  # Set up table for survivor output
  survivors <- vector("list", 100)
  pop_size <- vector("list", 100) # Right now this is hard-coded to count population size each year for 100 years. Will need to update this when I update the model run-time for dynamic populations (not just a single cohort that all die in under 100 years)
  
  # As long as more than 10 people are alive
  while (sum(!pop$dead) >= 10) {
    k <- k + 1  # Increment time
    pop <- age_pop(pop, k)
    Alive <- which(!pop$dead)
    
    # Baseline Siler mortality risk for this age (vector if the population is age-structured; scalar if the population is an age cohort)
age_based_risk <- compute_siler_risk(ages = pop[Alive,]$age, mortality_regime = mortality_regime) ###### Need to update this function. It relies on k, the time marker, which is only equivalent to agent age when the model generates a cohort and not an age-structured population.
    
    # Vectorized lesion formation and mortality across all living agents.
    # Each agent receives its own independent draw; case_when() dispatches
    # each agent through the correct mortality branch element-wise.
    if(!is.null(lesion_formation_rate)){
      pop <- form_lesions(pop, Alive,
                                 formation_window_opens,
                                 formation_window_closes,
                                 lesion_formation_rate)
    }
    
    
    pop <- apply_mortality(pop, Alive,
                                  age_based_risk,
                                  mortality_risk_type,
                                  relative_mortality_risk)
    
    # Update summary log for survivors. NOTE: Right now, the 'survivors' log isn't set up for age-structured populations. For now, if the model generates an age-structured population, it records population size annually instead. It doesn't track lesion frequency in the population as a whole. This is something to come back to in future, but for now, the survivors frequency table isn't particularly important for any of the questions we're applying this model to. 
    if(age_structured == FALSE){
    survivors[[k]] <- record_cohort_survivors(pop, k, lesion_formation_rate)
    }else{
    pop_size[[k]] <- data.frame(year = k+1, n = length(Alive))
    }
 }
  # Once 10 or fewer people are left alive — they all enter the cemetery
  pop <- finalize_cemetery(pop, k)
  
  # Apply deposition bias (if any)
  pop <- apply_deposition(pop,
                             deposition_model = 'cutoff',
                             deposition_param = deposition_param,
                             dx = 1)
  pop$in_sample <- pop$was_deposited
  
  # Apply preservation bias (if any)
  if (loss_strength != 'no_decay') {
    a_siler <- c(taphonomy_regime$a1, taphonomy_regime$b1,
                 taphonomy_regime$a2, taphonomy_regime$a3,
                 taphonomy_regime$b3)
    b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
    pop <- apply_preservation(pop,
                                 preservation_model = 'siler',
                                 preservation_param = b_siler,
                                 dx = 1)
  }
  
  # Apply age misestimation (if any)
  if (age_noise) pop <- apply_estimation_error(pop)
  
  # Remove internal columns before returning
  pop <- pop %>% dplyr::select(-"dead")
  
  # Model output
  # for a cohort model
  if(age_structured == FALSE){
  output <- list(individual_outcomes = pop, survivors = rbind(c("Age" = 1, "Alive" = pop0_size),do.call(rbind, survivors)))
  } else{
  # for an age-structured population
    output <- list(individual_outcomes = pop, population_size = rbind(c("Year" = 1, "n" = pop0_size), do.call(rbind, pop_size)))
  }
  
  return(output)
}
