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
#' @param model_lesions Logical. If true, then a column for lesion presence is initialized; no one has lesions atcurrent_time = 0. 
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
      pop0 <- pop0 %>%
        mutate(lesion = if_else(pop0$age %in% lesion_formation_window[1]:lesion_formation_window[2], 0, NA_real_)) %>% 
        relocate(lesion, .after = age) # change position of lesion column so it sits to the right of 'age'
      }
   
     return(pop0)
    }


#' Age up all living agents to the current timestep
#' @param pop Population data frame
#' @return Updated pop data frame
#' @keywords internal
age_pop <- function(pop) {
  pop$age <- pop$age + 1
  pop
}

#### @Saige: demohaz robust Siler parameters go here? Or we can just save the mortality regime parameter values as the robust ones. But for most users, I think we want to build in a transformation so they don't accidentally feed the wrong parameterization into the model and have it fail silently. 
#' Compute Siler hazard for a given age
#' @param ages Integer vector, if population is age-structured
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric hazard value
#' @keywords internal
compute_siler_risk <- function(ages, mortality_regime) { 
  if (any(mortality_regime < 0)) {
    stop('No parameters of the mortality regime can be negative')
  }
  if(mortality_regime$a3 < 1){ # These are traditional Siler parameters
    age_based_risk <- mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }
  if(mortality_regime$a3 > 1){ # These are robust Siler parameters
    mortality_regime <- demohaz_to_trad_siler_param(mortality_regime)
    age_based_risk <- mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }
  age_based_risk
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



#' Apply fertility to the living population for onecurrent_timestep
#'
#' Determines how many births occur thiscurrent_timestep and appends new agents to
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
#' @param tfr Numeric. Total Fertility Rate.
#' @param lesions Anything, or Null. If != Null, the lesions column is initialized with 0's. 
#' @param asfr Named numeric vector of age-specific fertility rates, as
#'   returned by compute_trapezoid_asfr(). The names must be character
#'   representations of integer ages.
#' @param dx Numeric.current_timestep size. Scales fertility rates proportionally.
#'   Default 1 (annualcurrent_timestep).
#' @return Updated pop data frame with new agents appended
#' @keywords internal
apply_fertility <- function(pop, tfr, asfr, dx = 1, model_lesions) {
  
  repro_ages      <- as.numeric(names(asfr))
  in_repro_window <- pop$age %in% repro_ages
  repro_ages_actual <- pop$age[in_repro_window]
  
  n_repro         <- sum(in_repro_window)
  n_female_approx <- round(n_repro / 2)
  
  if (n_female_approx == 0L) return(pop)
  
  female_age_sample <- sample(repro_ages_actual,
                              size    = n_female_approx,
                              replace = FALSE)
  
  asfr_values     <- asfr[as.character(female_age_sample)]
  expected_births <- sum(asfr_values * dx)
  n_births        <- rpois(1, lambda = expected_births)
  
  new_agents <- generate_births(pop, n_births, model_lesions)
  rbind(pop, new_agents)
}
  


#' Vectorized lesion formation across all living agents
#'
#' Replaces the per-agent form_lesion() loop with a single vectorized draw.
#' Each living agent gets its own uniform draw; conditions are evaluated
#' element-wise so individual-level heterogeneity is preserved.
#'
#' @param pop Population data frame
#' @param formation_window_opens Age at which lesions can start forming
#' @param formation_window_closes Age at which lesions stop forming
#' @param lesion_formation_rate Annual probability of lesion formation
#' @return Updated pop data frame
#' @keywords internal
form_lesions <- function(pop,
                         formation_window_opens,
                         formation_window_closes,
                         lesion_formation_rate) {
  stress <- runif(nrow(pop), 0, 1)
  
  in_window  <- pop$age >= formation_window_opens & pop$age <= formation_window_closes
  new_lesion <- as.integer(in_window & stress <= lesion_formation_rate)
  pop$lesion <- pmax(pop$lesion, new_lesion, na.rm = TRUE)
  
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
apply_mortality <- function(pop, 
                            age_based_risk,
                            mortality_risk_type,
                            relative_mortality_risk) {
  death_dice <- runif(nrow(pop), 0, 1)
  
  ages   <- pop$age
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
  
  pop$dead <- pop$dead | (death_dice < threshold)
  
  pop
}

#' Record survivor snapshot for the currentcurrent_timestep. This function is only built to handle a cohort, not a full age-structured population. *This could/should be changed*
#' @param pop Population data frame
#' @paramcurrent_time Currentcurrent_timestep
#' @param model_lesions Logical
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal
record_cohort_survivors <- function(pop, current_time, model_lesions) {
  n_alive  <- sum(!pop$dead & pop$age == current_time)
  n_lesion <- sum(!pop$dead & pop$lesion == 1 & pop$age == current_time)
  if(model_lesions){
    data.frame(Age = current_time, # in an age-cohort, time is equivalent to age. 
               Alive = n_alive,
               Lesion = n_lesion,
               Lesion_perc = ifelse(n_alive == 0, NA, round(n_lesion / n_alive * 100, 1)))
  }
  else{
    data.frame(Age = current_time,
               Alive = n_alive)
  }
}

#' Finalize the pop into a cemetery
#'
#' Ages remaining survivors one final step and marks all as dead.
#' @param pop Population data frame
#' @paramcurrent_time Last completed timestep
#' @return Data frame with all agents marked dead
#' @keywords internal
finalize_cemetery <- function(pop, decedents, current_time) {
 current_time <- current_time + 1
  if(nrow(pop) > 10){
    print("Too many agents are still alive -- don't kill them all yet!")
  }
  else{
    pop$age <- age_pop(pop) 
    pop$dead <- "dead"
    decedents[[current_time]] <- pop[,! colnames(pop) %in% c("dead") ]
    decedents[[current_time]]$year_died <- current_time 
  }
}


# --- Main simulation function (exported) ---

#' Simulate a cemetery using the Persephone ABM
#'
#' Runs the Persephone agent-based model: a birth pop ages throughcurrent_time,
#' facing annual risks of skeletal lesion formation and Siler-model mortality.
#' Individuals with lesions may experience modified mortality risk. The
#' simulation ends when fewer than 10 individuals remain alive.
#'
#' @param pop0_size Integer. Number of individuals in the starting pop.
#' @param dx Numeric. Size of thecurrent_timestep in which all other model actions are applied.
#' @param max_years. Integer. Number of years to run the simulation, given that the population doesn't crash before then. 
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
#' @param deposition_param
#' @param taphonomy_regime
#' @param loss_strength
#' @param age_noise
#' @return A list with two elements
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
                              formation_window_opens = 0, #### @ I'd like to change this to a c(0,0) style vector called 'lesion_formation_window'
                              formation_window_closes = 0,
                              mortality_risk_type = "proportional",
                              relative_mortality_risk = 1,
                              tfr, # @ haven't tried this with non-integer values yet.
                              mortality_regime,
                              pop_growth_rate = 0, # defaults to stationary population
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
  
current_time <- 0  # Initializecurrent_time counter
  
  # Generate starting cohort or age-structured population attime = 0.
  pop <- create_pop(pop0_size, age_structured = age_structured, 
                    r = pop_growth_rate, 
                    model_lesions = model_lesions,
                    lesion_formation_window = c(formation_window_opens, formation_window_closes),
                    mortality_regime = mortality_regime)
  
  # Calculate age-specific fertility rate, based on total fertility rate.
  if(age_structured == TRUE){
    asfr <- compute_trapezoid_asfr(ages = 0:100, tfr = tfr)
  }

  # Set up lists to track output
  survivors <- vector("list", max_years)
  pop_size <- vector("list", max_years) 
  decedents <- vector("list", max_years)
  
  
  # As long as more than 10 people are alive
  while (nrow(pop) > 10) {
    
    # break the loop if thecurrent_time has exceeded max_years
    if(current_time > max_years){
      break 
    }
    
   current_time <-current_time + 1
    pop <- age_pop(pop, current_time)
    
    age_based_risk <- compute_siler_risk(ages = pop$age, mortality_regime = mortality_regime)
    
    # Bones change.
    if (model_lesions) {
      pop <- form_lesions(pop,
                          formation_window_opens,
                          formation_window_closes,
                          lesion_formation_rate)
    }
    
    # The reaper comes. 
    pop <- apply_mortality(pop,
                           age_based_risk,
                           mortality_risk_type,
                           relative_mortality_risk)
    
    # Bring out yer dead. 
    decedents[[current_time]] <- pop[pop$dead, !colnames(pop) %in% c("dead")]
    if (nrow(decedents[[current_time]]) > 0) {
      decedents[[current_time]]$year_died <- current_time 
    }
    
    pop <- pop[!pop$dead, ]
    
    if (age_structured == TRUE) {
      
      # The stork visits. 
      pop <- apply_fertility(pop, tfr = tfr,
                             asfr = asfr, dx = dx, model_lesions = model_lesions)
    }
    
    if (age_structured == FALSE) {
      survivors[[current_time]] <- record_cohort_survivors(pop, current_time, lesion_formation_rate)
    } else {
      pop_size[[current_time]] <- data.frame(year = current_time + 1, n = nrow(pop))
    }
  }
  
  if(current_time <= max_years){
  # Once 10 or fewer people are left alive — they all enter the cemetery
  decedents[[k]] <- finalize_cemetery(pop, decedents, current_time )
  }
  
  # Apply deposition bias (if any)
  # decedents <- apply_deposition(decedents,
  #                            deposition_model = 'cutoff',
  #                            deposition_param = deposition_param,
  #                            dx = 1)
  # decedents$in_sample <- decedents$was_deposited
  # 
  # Apply preservation bias (if any)
  if (loss_strength != 'no_decay') {
    a_siler <- c(taphonomy_regime$a1, taphonomy_regime$b1,
                 taphonomy_regime$a2, taphonomy_regime$a3,
                 taphonomy_regime$b3)
    b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
    decedents <- apply_preservation(decedents,
                                 preservation_model = 'siler',
                                 preservation_param = b_siler,
                                 dx = 1)
  }
  
  # Apply age misestimation (if any)
  if (age_noise) decedents <- apply_estimation_error(decedents)
  
  
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
