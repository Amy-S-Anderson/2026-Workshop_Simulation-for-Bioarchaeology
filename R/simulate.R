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
                              lesion_related_hazard = 1, # This is only called if lesion modifies mortality hazard directly. I think this should be removed, since exposure itself is the thing that causes the hazard. The lesion is only an indicator of exposure to lesion-causing conditions. 

                              # Frailty arguments
                              gammafrailty_variance = NULL,
                              

                              # Exposure-lesion-hazard relationships
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
  
  # Check Input: return message if user chooses nonsensical combination of argument values.
  if (!is.null(lesion_formation_rate) && lesion_formation_window[2] == 0) {
    warning("lesion_formation_rate is set but formation window closes at 0: no lesions will form.")
  }
  if(is.numeric(tfr) && age_structured == FALSE){
    warning("fertility only works for an age-structured population. Either run a cohort with age_structured = FALSE and tfr = NULL, or run an age-structured population with age_structured = TRUE and tfr = a numeric value (0-12). Note that run time may be slow if tfr is large, especially if pop0_size and max_years are also on the high end for archaeological contexts.")
  }
  if(is.null(tfr) && age_structured == TRUE){
    warning("You are modeling an age-structured population but have not entered a fertility rate for the tfr argument. Please choose a numeric value for tfr, or change age_structured from TRUE to FALSE to model an age cohort.")
  }
  
  # Logical: Does this simulation include skeletal lesion formation?
  model_lesions <- !is.null(lesion_formation_rate) | !is.null(annual_exposure)
  # Logical: Does this simulation include heterogeneous frailty as a stable characteristic of agents, assigned at birth?

  # Hold population configuration parameters here for ease of reference in functions that follow. 
  pop_config <- list(
    model_lesions              = model_lesions,
    annual_exposure            = annual_exposure,
    lesion_formation_window.   = lesion_formation_window,
    frailty_variance           = gammafrailty_variance
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
    # If there is variance in assigned frailty values at birth, defrail the Siler values.
    defrail_siler(mortality_regime = params$mortality_regime,
                  frailty_variance  = pop_config$frailty_variance)
  } else {
    # Otherwise, use the Siler function directly to calculate age-based mortality hazards for ages 0 to 110. 
    ages <- 0:110
    data.frame(
      age = ages,
      mu0 = compute_siler_risk(ages, params$mortality_regime)
    )
  }
  # Add a column to the age-based hazards lookup table detailing the proportion of deaths at each age due to trauma.
  # This column incorporates life history differences in relative risk of dying from illness/natural causes vs. accidental/traumatic causes
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
  while (nrow(pop) > 0){
    # exit the loop when current_time is greater than max_years, or when the population drops to <=10 people. 
    if(current_time >= max_years) break
    # if <=10 people are alive, they will all die this time step. This truncation decision is based on the poor precision/accuracy of skeletal age-at-death estimates at old ages, and to prevent a stochastic model from producing an age outlier; bioarchaeologically we wouldn't be able to see Methuselah.
    force_death <- nrow(pop) <= 10  
    
    if(hazard_is_transient && "acquired_frailty" %in% names(pop)){ ### is this the right way to do a dual if statement? Maybe check your book on R code for if conditionals. 
      # re-set the value for acquired frailty back to baseline if hazard is transient. 
      pop$acquired_frailty <- NA_real_
    }
    
    # Sample agent exposure (to unspecified (stress?) event) first, so both form_lesions() and apply_mortality() can read it.
    # Determine who has experienced a 'stress' event in this time step. Flag them, and update their stress event count. 
    if (!is.null(annual_exposure) || !is.null(lesion_formation_rate)) {
      pop <- sample_exposure(pop, annual_exposure, lesion_formation_rate, lesion_formation_window) 
      # pop now has columns 'exposed_this_step' and 'n_stress_events'. 
      # n_stress_events already includes the current stress event in its total n. 
     # cat("year:", current_time, " n:", nrow(pop), " age range:", range(pop$age), "\n")
    }
    
   
    # If lesions are modeled in this simulation, then:
    if(model_lesions) {
      
    # Lesion formation timing is an open question that introduces an order-of-operations issue:
      # **Option 1: Lesions require a period of survival in order to form.** 
      # Apply mortality first, then form lesions among survivors. 
    if (lesion_requires_survival) {
      # Mortality first; lesion forms only for survivors
      updated_pop <- apply_mortality(pop,
                                     mu0_lookup = mu0_table,
                                     mortality_risk_type,
                                     force_death = force_death, # TRUE if 10 or fewer people are left alive in this time step. 
                                     current_time            = current_time,
                                     risk_factors            = list(
                                       frailty = NULL, # NULL here means that the actual value in the column 'frailty' will be used as a hazard multiplier. 
                                       acquired_frailty = NULL), 
                                     exposure_hazard_multiplier = exposure_hazard_multiplier) # a scalar.
      
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
      # Now apply lesion formation to survivors. 
      pop <- form_lesions(pop,
                            lesion_formation_window,
                            lesion_formation_rate,
                            annual_exposure)
      }
    
     else { # Lesion does not require a period of survival.
       
      # **Option 2:Lesions form immediately.** 
      # Form lesions.
      pop <- form_lesions(pop,
                          lesion_formation_window,
                          lesion_formation_rate,
                          annual_exposure)
      # Then apply mortality (below).
     }
    }
      # If lesions are not being modeled in this simulation, then skip lesion formation.
      # Jump straight to mortality.
      updated_pop <- apply_mortality(pop,
                                     mu0_lookup = mu0_table,
                                     mortality_risk_type,
                                     force_death = force_death, 
                                     current_time            = current_time,
                                     risk_factors            = list(
                                       frailty = NULL, 
                                       acquired_frailty = NULL), 
                                     exposure_hazard_multiplier = exposure_hazard_multiplier) 
                                     
      # and update the bookkeeping of the living and the dead. 
      decedents[[current_time]] <- updated_pop$decedents
      pop <- updated_pop$pop
    
    
    # If this is a dynamic population with a specified TFR for fertility, then apply fertility.
    # i.e., calculate new births and generate new agents.
    if (age_structured == TRUE & !is.null(tfr) && nrow(pop) > 0) {
      pop <- apply_fertility(pop, tfr = tfr, asfr = asfr, time = current_time, dx = dx, pop_config)
    }
      
    # Happy New Year!
    current_time <- current_time + 1
    
    # Happy Birthday!
    pop <- age_pop(pop)
    
    # If this is an age cohort, no fertility is in play. 
    # Record the time step/agent age (equivalent), and how many agents have survived this year.
    # Also record how many surviving agents have lesions, if lesions are being modeled in this simulation.
    if (age_structured == FALSE && nrow(pop) > 0) {
      survivors[[current_time]] <- record_cohort_survivors(pop, current_time, model_lesions)
    } else {
    # Or, if this is a dynamic population, record the year/time step, and the number of agents alive at this time step.
    # Do not record lesion counts, since lesion counts are sensitive to population age structure, and the number is not easily interpretable. 
      pop_size[[current_time]] <- data.frame(Time = current_time, n = nrow(pop))
    }
    if(yearly_updates){
      print(paste0("Year ", current_time))
      print(paste0("Pop_size = ", nrow(pop)))
    }
  }
  
  
  # Compile the Cemetery: Transform decedents from a list of data frames to a single data frame
  decedents <- do.call(rbind, decedents)

  #--------------------------------------------------------------------
  # Post-mortem Processes
  #--------------------------------------------------------------------

  # Apply deposition bias (if any)
  # This is a way to operationalize a common mortuary behavior; the youngest individuals are often buried in homes or places other than the cemetery where most others are buried.
  decedents <- apply_deposition(decedents,
                             deposition_model = 'cutoff',
                             deposition_param = deposition_param,
                             dx = 1)
  # create 'in_sample' variable. Populate with TRUE values for individuals who were deposited in the cemetery.
  decedents$in_sample <- decedents$was_deposited

  # Apply preservation bias (if any)
  # The youngest and oldest individuals are less well preserved, so we can model preservation bias using a Siler function, like mortality hazard.
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
  if (age_noise) {
    decedents <- apply_estimation_error(decedents)
  }

  
  # Organize model output
    if(model_lesions){
      starting_cohort <- c("Time" = 1, "Age" = 0, "n" = pop0_size, "Lesion" = 0, "Lesion_perc" = 0.0)
    } else{
      starting_cohort <- c("Time" = 1, "Age" = 0, "n" = pop0_size)
    }
  if(!is.null(tfr)){ 
    # If it is a dynamic population, record its size each year.
    annual_census = do.call(rbind, pop_size)
  } else{ 
    # If it is a single generation of agents, record the survivors in each year.
    annual_census = rbind(starting_cohort, do.call(rbind, survivors))
  }
  # Output:
  # 1. individual_outcomes, a simulated bioarchaeological data set 
  # 2. annual_census, the yearly record of number of living individuals
  output <- list(individual_outcomes = decedents, 
                 annual_census = annual_census)
  
  return(output)
}


