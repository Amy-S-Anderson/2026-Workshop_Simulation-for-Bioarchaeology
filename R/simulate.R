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
                                     max_age = 100) {
  ages <- 0:max_age
  
  sad <- compute_stable_age_distribution(ages, mortality_regime, r)
  
  # Sample each agent's starting age from the stable age distribution
  starting_ages <- sample(ages, size = pop0_size, replace = TRUE, prob = sad)
  
   data.frame(
    agent_id     = 1:pop0_size,
    age          = as.numeric(starting_ages),
    dead         = FALSE,
    was_deposited = FALSE,
    in_sample    = TRUE
  )
}

#' @param pop_size Number of agents in the starting pop
#' @param age_structured Logical: is this an age-structured population? (if not, it is a single age-cohort)
#' @param model_lesions Logical. If true, then a column for lesion presence is initialized; no one has lesions at time = 0. 
#' @param lesion_formation_window Vector length 2: c(age at which window opens, age at which window closes) 
#' @param r Numeric, the population growth rate
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return A data frame with columns: agent_id, age, lesion, dead, in_sample
#' @keywords internal
create_pop <- function(pop0_size, age_structured, 
                       model_lesions, 
                       lesion_formation_window,
                       r = 0, mortality_regime) {
  if(age_structured == TRUE){
    pop0 <- create_pop_stable_age(pop0_size,
                                  mortality_regime,
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
  }
    if(model_lesions){
      pop0$lesion <- if_else(pop0$age %in% lesion_formation_window[1]:lesion_formation_window[2], 0, NA)
      }
    pop0 <- pop0 %>% relocate(lesion, .after = age) # change position of lesion column so it sits to the right of 'age'
   
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
    mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }


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
compute_trapezoid_asfr <- function(ages, tfr,
                                   age_start = 15, age_peak_start = 25,
                                   age_peak_end = 35, age_end = 45) {
  flat_width      <- age_peak_end   - age_peak_start
  ramp_up_width   <- age_peak_start - age_start
  ramp_down_width <- age_end        - age_peak_end
  area_per_peak   <- flat_width + ramp_up_width / 2 + ramp_down_width / 2
  peak_asfr       <- tfr / area_per_peak
  
  asfr <- dplyr::case_when(
    ages <  age_start      ~ 0,
    ages <  age_peak_start ~ peak_asfr * (ages - age_start)  / ramp_up_width,
    ages <= age_peak_end   ~ peak_asfr,
    ages <  age_end        ~ peak_asfr * (age_end - ages)    / ramp_down_width,
    TRUE               ~ 0
  )
  
  setNames(asfr, ages)
}


#' Generate new agent rows to append to the pop
#'
#' Creates n_births new agents with age 0 and no lesion, with agent_ids
#' continuing sequentially from the current maximum in the population. All other
#' columns are initialized to match the structure of create_pop().
#'
#' @param pop Population data frame (used to determine next agent_id)
#' @param n_births Integer. Number of new agents to create.
#' @param model_lesions Logical. If FALSE, no lesion column is initalized. 
#' @return Data frame with n_births rows ready to rbind() onto the population
#' @keywords internal
generate_births <- function(pop, n_births, model_lesions) {
  if (n_births == 0L) return(pop[0L, ])  # empty frame with correct columns
  
  next_id <- max(pop$agent_id) + 1L
  babies <- data.frame(agent_id = seq(next_id, next_id + n_births - 1L),
                     age = 0,
                     dead = FALSE,
                     was_deposited = FALSE,
                     in_sample = TRUE
  )
  if(model_lesions){
    babies <- babies %>%
      mutate(lesion = 0L) %>%
      relocate(lesion, .after = age)
  }
  return (babies)
}



#' Apply fertility to the living population for one timestep
#'
#' Determines how many births occur this timestep and appends new agents to
#' the pop. The number of births is drawn from a Poisson distribution whose
#' rate equals the sum of age-specific fertility rates across all reproductive-
#' age women in the living population.
#'
#' Women are approximated as half of all living agents aged between age_start
#' and age_end (inclusive). This assumption can be replaced once sex is added
#' as a column to the pop data frame, by filtering on pop$sex == "F"
#' instead.
#'
#' Birth counts are drawn as Poisson rather than Bernoulli because multiple
#' births can occur at a given age in a population of any size, and the Poisson
#' approximation to the sum of many independent Bernoulli trials is exact in
#' the limit of large population and small per-individual rate.
#'
#' @param pop Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param tfr Numeric. Total Fertility Rate.
#' @param lesions Anything, or Null. If != Null, the lesions column is initialized with 0's. 
#' @param asfr Named numeric vector of age-specific fertility rates, as
#'   returned by compute_trapezoid_asfr(). The names must be character
#'   representations of integer ages.
#' @param dx Numeric. Timestep size. Scales fertility rates proportionally.
#'   Default 1 (annual timestep).
#' @return Updated pop data frame with new agents appended
#' @keywords internal
apply_fertility <- function(pop, Alive, tfr, asfr, dx = 1, model_lesions) {
  
  # Ages of all living agents
  ages_alive <- pop$age[Alive]
  
  # Identify reproductive-age agents and approximate females as half of them
  repro_ages    <- as.numeric(names(asfr))
  in_repro_window <- ages_alive %in% repro_ages
  
  repro_agents  <- Alive[in_repro_window]
  repro_ages_actual <- pop$age[repro_agents]
  
  # Approximate females as half the reproductive-age living population.
  # When sex becomes a pop column, replace this block with:
  #   female_agents <- repro_agents[pop$sex[repro_agents] == "F"]
  #   female_ages   <- pop$age[female_agents]
  n_repro       <- length(repro_agents)
  n_female_approx <- round(n_repro / 2)
  
  if (n_female_approx == 0L) return(pop)  # no reproductive-age women
  
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
  
  # Append new agents to the pop and return
  new_agents <- generate_births(pop, n_births, model_lesions)
  pop <- rbind(pop, new_agents)
  
  pop
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
  pop$lesion[Alive] <- pmax(pop$lesion[Alive], new_lesion, na.rm = TRUE)
  
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
  has_lesions <- "lesion" %in% names(pop) # does the lesion column exist in this simulation?
  
  # Compute the effective death threshold for each agent.
  # case_when() evaluates conditions row-wise (element-wise on vectors), so
  # each agent follows exactly one branch — logically equivalent to the
  # original nested ifelse() but readable and easy to extend.
  if(has_lesions){
  threshold <- dplyr::case_when(
    pop$lesion == 0 ~
      age_based_risk,
    
    pop$lesion == 1 & mortality_risk_type == "proportional" ~
      age_based_risk * relative_mortality_risk,
    
    pop$lesion == 1 & mortality_risk_type == "time_decreasing" ~
      age_based_risk * relative_mortality_risk / ((ages / 10) + relative_mortality_risk),
    
    pop$lesion == 1 & mortality_risk_type == "time_increasing" ~
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
#' @param model_lesions Logical
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal
record_cohort_survivors <- function(pop, k, model_lesions) {
  n_alive  <- sum(!pop$dead & pop$age == k)
  n_lesion <- sum(!pop$dead & pop$lesion == 1 & pop$age == k)
  if(model_lesions){
    data.frame(Age = k,
               Alive = n_alive,
               Lesion = n_lesion,
               Lesion_perc = ifelse(n_alive == 0, NA, round(n_lesion / n_alive * 100, 1)))
  }
  else{
    data.frame(Age = k,
               Alive = n_alive)
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
#' @param tfr Numeric. Total fertility rate. 
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
                              dx = 1, 
                              max_years = 100,
                              lesion_formation_rate = NULL,
                              formation_window_opens = 0,
                              formation_window_closes = 0,
                              mortality_risk_type = "proportional",
                              relative_mortality_risk = 1,
                              tfr,
                              mortality_regime,
                              pop_growth_rate = 0, # default to stationary population
                              age_structured = TRUE,
                              deposition_param = 0,
                              taphonomy_regime,
                              loss_strength = 'no_decay',
                              age_noise = FALSE) {
  
  # check: return message if user chooses nonsensical combination of argument values.
  if (!is.null(lesion_formation_rate) && formation_window_closes == 0) {
    warning("lesion_formation_rate is set but formation_window_closes is 0: no lesions will form.")
  }
  
  # Logical: Does this simulation include skeletal lesion formation?
  model_lesions <- !is.null(lesion_formation_rate)
  
  k <- 0  # Initialize time counter
  
  # Generate starting cohort or age-structured population at k = 0.
  pop <- create_pop(pop0_size, age_structured = age_structured, 
                    r = pop_growth_rate, 
                    model_lesions = model_lesions,
                    lesion_formation_window = c(formation_window_opens, formation_window_closes),
                    mortality_regime = mortality_regime)
  
  # Calculate age-specific fertility rate, based on total fertility rate.
  if(age_structured == TRUE){
    asfr <- compute_trapezoid_asfr(ages = 0:100, tfr = tfr)
  }

  # Set up table for survivor output
  survivors <- vector("list", 100)
  pop_size <- vector("list", 100) # Right now this is hard-coded to count population size each year for 100 years. Will need to update this when I update the model run-time for dynamic populations (not just a single cohort that all die in under 100 years)
  decedents <- vector("list", max_years)
  
  
  # As long as more than 10 people are alive
  while (sum(!pop$dead) >= 10) {
    k <- k + 1  # Increment time
    pop <- age_pop(pop, k)
    Alive <- which(!pop$dead)
    
    # Baseline age-dependent Siler mortality risk for all individuals (vector if the population is age-structured; scalar if the population is an age cohort)
age_based_risk <- compute_siler_risk(ages = pop[Alive,]$age, mortality_regime = mortality_regime) ###### Need to update this function. It relies on k, the time marker, which is only equivalent to agent age when the model generates a cohort and not an age-structured population.
    
    # Vectorized lesion formation and mortality across all living agents.
    # Each agent receives its own independent draw; case_when() dispatches
    # each agent through the correct mortality branch element-wise.
    if(model_lesions){
      
      # Bones change.
      pop <- form_lesions(pop, Alive,
                                 formation_window_opens,
                                 formation_window_closes,
                                 lesion_formation_rate)
    }
    
    # The Reaper comes.
    pop <- apply_mortality(pop, Alive,
                                  age_based_risk,
                                  mortality_risk_type,
                                  relative_mortality_risk)
    
    # Bring out yer dead. 
    decedents[[k]]        <- pop[pop$dead, ]
    if(nrow(decedents[[k]]) > 0){ # if anyone died
      decedents[[k]]$year_died <- k # record the year of their death
    }

    # Separate the dead from the living. 
    pop                  <- pop[!pop$dead, ]
    
    if(age_structured == TRUE){
      # The stork visits.
    pop                  <- apply_fertility(pop, Alive, tfr = tfr, 
                                            asfr = asfr, dx = dx, model_lesions = model_lesions)  
    }
    # NOTE: Every call to apply_fertility() does rbind(pop, new_agents), which copies the entire population data frame. For a long-running simulation with high fertility this compounds quickly. The same list-then-bind pattern you've already applied to survivors would help here: accumulate new birth batches in a list and bind them to pop every N timesteps or at the end of the simulation. That said, this is a lower priority optimization than the others — worth flagging for when you're tuning for large populations.
    
    
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
  output <- list(individual_outcomes = do.call(rbind, decedents), survivors = rbind(c("Age" = 1, "Alive" = pop0_size),do.call(rbind, survivors)))
  } else{
  # for an age-structured population
    output <- list(individual_outcomes = do.call(rbind, decedents), population_size = rbind(c("Year" = 1, "n" = pop0_size), do.call(rbind, pop_size)))
  }
  
  return(output)
}
