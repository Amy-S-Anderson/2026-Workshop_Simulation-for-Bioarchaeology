# --------------------------------------------------------------------------------------------------------------
# PERSEPHONE: Agent-Based Model for Bioarchaeology
# --------------------------------------------------------------------------------------------------------------
#
# This file contains the core simulation engine and its helper functions.
# Individual agents:
#   - Are born as a starting pop of a specified size
#   - Experience a risk of dying in each time step of the model
#   - Age-specific mortality risks follow a Siler function


# Depending on other option arguments, agents may:
#   - Face a stable annual risk of forming a skeletal lesion within a specified age range
#   - Experience different age-specific mortality depending on lesion status
#   - Experience different age-specific risks of traumatic death, and age/lesion-specific risks of illness-related death. 
#   - Be deposited (or not) in a cemetery. 
#   - Be preserved (or not) in the cemetery until time of recovery.
#   - Have their age-at-death estimated from skeletal features (which adds noise and bias to actual ages)






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
#' @param trauma_regime Data frame with key ages and the proportion of deaths due to traumatic injury at these ages. Defaults to NULL (no competing hazards in the model, just all-cause mortality)
#' @param pop_growth_rate Numeric. The population growth rate
#' @param age_structured Logical. Is the starting population age-structured, or an age cohort?
#' @param deposition_param Numeric. the minimum age for including agents in the cemetery
#' @param taphonomy_regime A named data frame storing Siler values. Taphonomic loss, like mortality hazard, is greatest for infants and elders, so we chose to model it using a Siler function.  
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
                              pop0_growth_rate = 0, # defaults to stationary population
                              tfr, # total fertility rate, on average, per woman
                              mortality_regime, # A named vector of Siler mortality hazard parameter values
                              trauma_regime = NULL, # A named data frame of key ages and the proportion of deaths due to traumatic injuries at these ages (gets at life history changes in competing hazards for deaths from injury vs. illness)
                              
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
  if(is.numeric(tfr) && age_structured == FALSE){
    warning("fertility only works for an age-structured population. Either run a cohort with age_structured = FALSE and tfr = NULL, or run an age-structured population with age_structured = TRUE and tfr = a numeric value (0-12). Note that run time may be slow if tfr is large, especially if pop0_size and max_years are also on the high end for archaeological contexts.")
  }
  
  # Logical: Does this simulation include skeletal lesion formation?
  model_lesions <- !is.null(lesion_formation_rate) | !is.null(annual_exposure)
  model_frailty <- !is.null(gammafrailty_variance)
  
  # Hold population configuration parameters here for ease of reference in functions that follow. 
  pop_config <- list(
    model_lesions              = model_lesions,
    annual_exposure            = annual_exposure,
    model_frailty              = model_frailty,
    frailty_variance           = gammafrailty_variance,
    lesion_formation_window    = lesion_formation_window,
    exposure_causes_hazard     = exposure_causes_hazard,
    hazard_is_transient        = hazard_is_transient,
    lesion_requires_survival   = lesion_requires_survival,
    exposure_hazard_multiplier = exposure_hazard_multiplier
  )
  
  #--------------------------------------------------------------------
  # Set up the Starting Population
  #--------------------------------------------------------------------
  # Generate starting cohort or age-structured population at time = 0.
  pop <- create_pop(pop0_size, age_structured = age_structured, 
                    r = pop0_growth_rate, 
                    mortality_regime = mortality_regime,
                    pop_config = pop_config)

  # Calculate age-specific fertility rate, based on total fertility rate.
  if(age_structured == TRUE & is.numeric(tfr)){
    asfr <- compute_trapezoid_asfr(ages = 0:100, tfr = tfr)
  }

  # Calculate age-specific mortality hazards across individual frailty values by 'defrailing' the Siler function, which gives average mortality hazard at each age. 
  mu0_table <- if (!is.null(pop_config$frailty_variance) && pop_config$frailty_variance > 0) {
    defrail_siler(mortality_regime = params$mortality_regime,
                  frailty_variance  = pop_config$frailty_variance)
  } else {
    ages <- 0:110
    data.frame(
      age = ages,
      mu0 = compute_siler_risk(ages, params$mortality_regime)
    )
  }
  # Add a column to the lookup table detailing the proportion of deaths at each age due to trauma (rather than illness)
  mu0_table <- build_p_trauma_lookup(mu0_lookup = mu0_table, control_points = trauma_regime)
  
  
  #--------------------------------------------------------------------
  # Get ready to run the clock forward and track the population
  #--------------------------------------------------------------------
  # Set up lists to track output
  survivors <- vector("list", max_years)
  pop_size <- vector("list", max_years) 
  decedents <- vector("list", max_years)
  
  
  current_time <- 1  # Initialize current_time counter
  
  #--------------------------------------------------------------------
  # Time begins!
  #--------------------------------------------------------------------
  # As long as people are alive
  while (nrow(pop) > 0) {
    # exit the loop when current_time is greater than max_years, or when the population drops to <=10 people. 
    force_death <- nrow(pop) <= 10  | (!is.null(tfr) && current_time > max_years) # max_years should only be checked  in a dynamic population, where TFR is specified and new agents are born into the model each year.
    
    
    # Sample exposure first, so both form_lesions() and apply_mortality() can read it.
    # Determine who has experienced a 'stress' event in this time step. Flag them, and update their stress event count. 
    if (!is.null(annual_exposure)) {
      pop <- sample_exposure(pop, annual_exposure) # pop now has columns 'exposed_this_step' and 'n_stress_events'
    }
    
    # If lesions are modeled in this simulation, then:
    if(model_lesions) {
    # Lesion formation timing is an open question that introduces an order-of-operations issue:
    # Option 1: Lesions require a period of survival in order to form. Apply mortality first, then form lesions among survivors. 
    if (lesion_requires_survival) {
      # Mortality first; lesion forms only for survivors
      updated_pop <- apply_mortality(pop,
                                     mu0_lookup = mu0_table,
                                     mortality_risk_type,
                                     force_death = force_death, # TRUE if 10 or fewer people are left alive in this time step. 
                                     current_time            = current_time,
                                     risk_factors            = list(
                                       frailty = NULL,
                                       acquired_frailty = NULL,
                                       lesion = lesion_related_hazard),
                                     lesion_requires_survival   = lesion_requires_survival,
                                     exposure_causes_hazard     = exposure_causes_hazard,
                                     hazard_is_transient        = hazard_is_transient,
                                     exposure_hazard_multiplier = exposure_hazard_multiplier)
      
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
      pop <- form_lesions(pop,
                            lesion_formation_window,
                            lesion_formation_rate,
                            annual_exposure,
                            exposure_causes_hazard     = exposure_causes_hazard,
                            hazard_is_transient        = hazard_is_transient,
                            exposure_hazard_multiplier = exposure_hazard_multiplier)
      
      pop$exposed_this_step <- NULL  # <-- re-set the values for immediate hazards here, after form_lesions() is done
      pop$transient_hazard  <- NULL  
    } else {
      # Option 2: Lesions don't require a period of survival after the lesion-causing event. Lesions form immediately, then mortality is calculated. 
        pop <- form_lesions(pop,
                            lesion_formation_window,
                            lesion_formation_rate,
                            annual_exposure,
                            exposure_causes_hazard     = exposure_causes_hazard,
                            hazard_is_transient        = hazard_is_transient,
                            exposure_hazard_multiplier = exposure_hazard_multiplier)
    }
  }
      
      updated_pop <- apply_mortality(pop,
                                     mu0_lookup = mu0_table,
                                     mortality_risk_type,
                                     force_death = force_death, # TRUE if 10 or fewer people are left alive in this time step. 
                                     current_time = current_time,
                                     risk_factors = list(frailty          = NULL,
                                                         acquired_frailty = NULL,
                                                         lesion           = lesion_related_hazard))
      
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
    }
    
    if (age_structured == TRUE & !is.null(tfr) && nrow(pop) > 0) {
      pop <- apply_fertility(pop, tfr = tfr, asfr = asfr, time = current_time, dx = dx, pop_config)
    }
    
    current_time <- current_time + 1
    pop <- age_pop(pop)
   # cat("Time:", current_time, "| nrow:", nrow(pop), "| age col exists:", "age" %in% names(pop), "| unique ages:", length(unique(pop$age)), "\n")
    
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
  # If this is a single-generation model (tfr = NULL, no agents born into the simulation), then When the population drops to <= 10 people, they all enter the cemetery.
  if(nrow(pop) > 0){
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
    # update the 'in_sample' variable to reflect individuals lost to taphonomic degradation
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
  # output = individual_outcomes, essentially a simulated bioarchaeological data set, and 
  # annual_census, the yearly record of number of living individuals (and, if a cohort and not an age-structured population and lesions are modeled in this simulation, also the number and % of individuals each year with a skeletal lesion)
  output <- list(individual_outcomes = decedents, 
                 annual_census = annual_census)
  
  return(output)
}


