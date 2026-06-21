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






#' Vectorized lesion formation across all living agents
#'
#' Replaces the per-agent form_lesion() loop with a single vectorized draw.
#' Each living agent gets its own uniform draw; conditions are evaluated
#' element-wise so individual-level heterogeneity is preserved.
#'
#' @param pop Population data frame
#' @param lesion_formation_window numeric vector of length = 2. c(Age at which lesions can start forming, Age at which lesions stop forming)
#' @param lesion_formation_rate A probability
#' @param annual_exposure A proportion. If lesion_formation_rate has a value then annual_exposure should be set to NULL, and vice versa. 
#' @param exposure_causes_hazard Logical. Does exposure to lesion-causing events affect mortality hazard?
#' @param hazard_is_transient Logical. If true, exposure-mediated change in mortality hazard last only for the year of exposure.
#' @param exposure_hazard_multiplier Numeric. The amount to multiply an agent's mortality hazard if exposure_causes_hazard == TRUE. 
#' @return Updated pop data frame
#' @keywords internal
#' Note: These aren't quite parallel options. if(lesion_formation_rate), then acquired frailty values aren't updated. Right now frailty is only integrated with deterministic lesion formation, not probabilistic lesion formation. 
form_lesions <- function(pop,
                         lesion_formation_window,
                         lesion_formation_rate = NULL,
                         annual_exposure = NULL,
                         exposure_causes_hazard = FALSE,
                         hazard_is_transient = FALSE,
                         exposure_hazard_multiplier = 1) {
  
  in_window <- pop$age >= lesion_formation_window[1] &
    pop$age <= lesion_formation_window[2]
  
  n_specified <- sum(!is.null(lesion_formation_rate), !is.null(annual_exposure))
  if (n_specified == 0 || n_specified == 2) {
    stop("Exactly one of lesion_formation_rate or annual_exposure must be specified.")
  }
  
  if (!is.null(lesion_formation_rate)) {
    # Original framework: probabilistic, no exposure concept
    stress <- runif(nrow(pop), 0, 1)
    exposed <- stress <= lesion_formation_rate
    pop$lesion <- pmax(pop$lesion, as.integer(in_window & exposed), na.rm = TRUE)
    return(pop)
  }
  
  # Annual exposure framework: reads exposed_this_step written by sample_exposure()
  exposed <- pop$exposed_this_step
  
  # Form lesions immediately (lesion_requires_survival = FALSE cases)
  pop$lesion <- pmax(pop$lesion, as.integer(in_window & exposed), na.rm = TRUE)
  
  # Apply hazard effects if applicable
  if (exposure_causes_hazard) {
    if (hazard_is_transient) {
      # Transient: write to transient_hazard, read by apply_mortality() this step only
      pop$transient_hazard <- ifelse(exposed, exposure_hazard_multiplier, 1)
    } else {
      # Permanent: accumulate into acquired_frailty
      if ("acquired_frailty" %in% names(pop)) {
        pop$acquired_frailty[is.na(pop$acquired_frailty)] <- 0
        pop$acquired_frailty <- ifelse(exposed,
                                       pop$acquired_frailty + exposure_hazard_multiplier,
                                       pop$acquired_frailty)
      }
    }
  }
  
  pop
}



 

#' Record survivor snapshot for the current timestep. 
#' @param pop Population data frame
#' @param current_time Current timestep
#' @param model_lesions Logical
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal
record_cohort_survivors <- function(pop, current_time, model_lesions) {
    survivors <- data.frame(Time = current_time,
                            Age = unique(pop$age), 
                            Alive = nrow(pop))
    if(model_lesions){
      survivors$Lesion = sum(pop$lesion, na.rm = TRUE)
      survivors$Lesion_perc = round((sum(pop$lesion, na.rm = TRUE) / nrow(pop)) * 100, 1)
    }
    
  survivors
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
#' @param lesion_formation_window  numeric vector of length = 2. c(Age at which lesions can start forming, Age at which lesions stop forming)
#' @param mortality_risk_type Character. How lesions modify mortality:
#'   "proportional", "time_decreasing", or "time_increasing".
#' @param lesion_related_hazard Numeric. Mortality multiplier for individuals
#'   with lesions. 1 = no effect, 2 = double risk. Default 1.
#' @param tfr Numeric. Total fertility rate. 
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#' @param pop_growth_rate Numeric. The population growth rate
#' @param age_structured Logical. Is the starting population age-structured, or an age cohort?
#' @param deposition_param Numeric. the minimum age for including agents in the cemetery
#' @param taphonomy_regime Not sure*** 
#' @param loss_strength Character. describes age-dependent preservation bias (defaults to 'no_decay')
#' @param age_noise Logical. If TRUE, then age estimation error is added to estimated age-at-death. 
#' @param yearly_updates Logical. If TRUE, then the model prints current year and population size at the end of each time step during the model run. Useful for assessing simulation progress, especially when timeframe, TFR, or population size are large.
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
#'   lesion_formation_window = c(0,5),
#'   mortality_regime = CoaleDemenyWestF5
#' )
#'
#' @export
Simulate_Cemetery <- function(# Time arguments
                              dx = 1, # size of time step in model time
                              max_years = 100, # max time this model will run, if the population doesn't crash first.
                              
                              # Demography arguments
                              pop0_size,
                              age_structured = TRUE, # if FALSE, this is a cohort model (proxy for stationary population)
                              pop_growth_rate = 0, # defaults to stationary population
                              tfr, # total fertility rate, on average, per woman
                              mortality_regime, # A named vector of Siler mortality hazard parameter values
                              
                              # Skeletal lesion arguments
                              lesion_formation_rate = NULL,
                              annual_exposure = NULL,
                              lesion_formation_window = c(0,0), 
                              mortality_risk_type = "proportional",
                              lesion_related_hazard = 1, # This is only called if lesion modifies mortality hazard directly

                              # Frailty arguments
                              gammafrailty_variance = NULL,

                              # Exposure-lesion-hazard relationships
                              exposure_causes_hazard   = FALSE,  # does exposure to lesion-causing events modify mortality?
                              hazard_is_transient      = FALSE,  # TRUE = year of exposure only, FALSE = permanent
                              lesion_requires_survival = FALSE,  # TRUE = agent must survive exposure year
                              exposure_hazard_multiplier = 1,    # the 'stress value' or acquired frailty that results from exposure to lesion-causing conditions. 
                              
                              # Post-mortem process arguments
                              deposition_param = 0,
                              taphonomy_regime = NULL,
                              loss_strength = 'no_decay',
                              age_noise = FALSE,
                              
                              # real-time update messages arguments
                              yearly_updates = FALSE) {
  
  # check: return message if user chooses nonsensical combination of argument values.
  if (!is.null(lesion_formation_rate) && lesion_formation_window[2] == 0) {
    warning("lesion_formation_rate is set but formation window closes at 0: no lesions will form.")
  }
  
  # Logical: Does this simulation include skeletal lesion formation?
  model_lesions <- !is.null(lesion_formation_rate) | !is.null(annual_exposure)
  model_frailty <- !is.null(gammafrailty_variance)
  
  # Hold population configuration parameters here for ease of reference in functions that follow. 
  pop_config <- list(
    model_lesions              = model_lesions,
    annual_exposure            = annual_exposure,
    model_frailty              = model_frailty,
    gammafrailty_variance         = gammafrailty_variance,
    lesion_formation_window    = lesion_formation_window,
    exposure_causes_hazard     = exposure_causes_hazard,
    hazard_is_transient        = hazard_is_transient,
    lesion_requires_survival   = lesion_requires_survival,
    exposure_hazard_multiplier = exposure_hazard_multiplier
  )
  

  # Generate starting cohort or age-structured population at time = 0.
  pop <- create_pop(pop0_size, age_structured = age_structured, 
                    r = pop_growth_rate, 
                    mortality_regime = mortality_regime,
                    pop_config = pop_config)

  # Calculate age-specific fertility rate, based on total fertility rate.
  if(age_structured == TRUE){
    asfr <- compute_trapezoid_asfr(ages = 0:100, tfr = tfr)
  }

  # Calculate age-specific mortality hazards across individual frailty values by 'defrailing' the Siler function, which gives average mortality hazard at each age. 
  mu0_table <- defrail_siler(mortality_regime, frailty_variance)
  
  # Set up lists to track output
  survivors <- vector("list", max_years)
  pop_size <- vector("list", max_years) 
  decedents <- vector("list", max_years)
  
  
  current_time <- 1  # Initialize current_time counter
  
  # As long as more than 10 people are alive
  while (nrow(pop) > 10) {
    
    # exit the loop when current_time is greater than max_years. This is only relevant when modeling a dynamic population, i.e., one where agents are born into the model world rather than just created by create_pop().
    if (!is.null(tfr) && current_time > max_years) break
    
    # Sample exposure first, so both form_lesions() and apply_mortality() can read it
    if (!is.null(annual_exposure)) {
      pop <- sample_exposure(pop, annual_exposure)
    }

    if (lesion_requires_survival) {
      # Mortality check first; lesion forms only for survivors
      updated_pop <- apply_mortality(pop,
                                     schedule = mortality_regime,
                                     mortality_risk_type,
                                     current_time            = current_time,
                                     risk_factors            = list(frailty          = NULL,
                                                                    acquired_frailty = NULL,
                                                                    lesion           = lesion_related_hazard),
                                     lesion_requires_survival   = lesion_requires_survival,
                                     exposure_causes_hazard     = exposure_causes_hazard,
                                     hazard_is_transient        = hazard_is_transient,
                                     exposure_hazard_multiplier = exposure_hazard_multiplier)
      
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
      
      if (model_lesions) {
        pop <- form_lesions(pop,
                            lesion_formation_window,
                            lesion_formation_rate,
                            annual_exposure,
                            exposure_causes_hazard     = exposure_causes_hazard,
                            hazard_is_transient        = hazard_is_transient,
                            exposure_hazard_multiplier = exposure_hazard_multiplier)
      }
      pop$exposed_this_step <- NULL  # <-- clean up here, after form_lesions() is done
      pop$transient_hazard  <- NULL  
    } else {
      # Lesion forms first, then mortality check
      if (model_lesions) {
        pop <- form_lesions(pop,
                            lesion_formation_window,
                            lesion_formation_rate,
                            annual_exposure,
                            exposure_causes_hazard     = exposure_causes_hazard,
                            hazard_is_transient        = hazard_is_transient,
                            exposure_hazard_multiplier = exposure_hazard_multiplier)
      }
      
      updated_pop <- apply_mortality(pop,
                                     schedule = mortality_regime,
                                     mortality_risk_type,
                                     current_time = current_time,
                                     risk_factors = list(frailty          = NULL,
                                                         acquired_frailty = NULL,
                                                         lesion           = lesion_related_hazard))
      
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
    }
    
    if (age_structured == TRUE) {
      pop <- apply_fertility(pop, tfr = tfr, asfr = asfr, time = current_time, dx = dx, pop_config)
    }
    
    current_time <- current_time + 1
    pop <- age_pop(pop)
    
    if (age_structured == FALSE) {
      survivors[[current_time]] <- record_cohort_survivors(pop, current_time, model_lesions)
    } else {
      pop_size[[current_time]] <- data.frame(Time = current_time, n = nrow(pop))
    }
    if(yearly_updates){
      print(paste0("Year ", current_time))
      print(paste0("Pop_size = ", nrow(pop)))
    }
  }
  
  # If this is a single-generation model (no Total Fertility Rate specified for agents)
  if(is.null(tfr)){
  # Once 10 or fewer people are left alive — they all enter the cemetery
  # This truncation decision is based on the poor precision/accuracy of skeletal age-at-death estimates at old ages, and to prevent a stochastic model from producing an age outlier; bioarchaeologically we wouldn't be able to see Methuselah. 
    pop <- age_pop(pop) 
    decedents[[current_time]] <- pop
    decedents[[current_time]]$year_died <- current_time 
  }
 
  # Transform decedents from a list of data frames to a single data frame
  decedents <- do.call(rbind, decedents) 
  
  # Add post-mortem processes to all decedents:
  # Apply deposition bias (if any)
  decedents <- apply_deposition(decedents,
                             deposition_model = 'cutoff',
                             deposition_param = deposition_param,
                             dx = 1)
  
  decedents$in_sample <- decedents$was_deposited

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
    if(model_lesions){
      starting_cohort <- c("Time" = 1, "Age" = 0, "n" = pop0_size, "Lesion" = 0, "Lesion_perc" = 0.0)
    } else{
      starting_cohort <- c("Time" = 1, "Age" = 0, "n" = pop0_size)
    }
  if(!is.null(tfr)){ # If it is a dynamic population, record its size each year
    annual_census = do.call(rbind, pop_size)
  } else{ # If it is a single generation of agents, record the survivors in each year
    annual_census = rbind(starting_cohort, do.call(rbind, survivors))
  }
  output <- list(individual_outcomes = decedents, 
                 annual_census = annual_census)
  
  return(output)
}


### Note: columns for 'was_deposited' and 'in_sample' need to be created at the end of the ABM function. They are no longer initialized in the starting population. 