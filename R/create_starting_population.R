#' @title Create a starting population or age-cohort of agents
#'
#' @description Generates a data frame of n agents with columns to hold specified agent states (age, lesion status, etc.), depending on the arguments passed to the function. 
#'
#' 
#' 
#' @param pop_size Number of agents in the starting pop
#' @param age_structured Logical: is this an age-structured population? (if not, it is a single age-cohort)
#' @param model_lesions Logical. If true, then a column for lesion presence is initialized; no one has lesions at current_time = 0. 
#' @param lesion_formation_window Vector length 2: c(age at which window opens, age at which window closes) 
#' @param gammafrailty_shape The alpha/shape parameter in a gamma distribution
#' @param gammafrailty_scale The sigma/scale parameter in a gamma distribution. Together these two parameters describe the frailty distribution at birth in the starting population.
#' @param r Numeric, the population growth rate
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @param pop_config List of population parameters
#' @return A data frame with columns: agent_id, age, lesion, dead, in_sample
#' @keywords internal
#' 
#' @export
create_pop <- function(pop0_size, age_structured, 
                       pop_config, # list of optional traits to initialize (lesions, frailty values)
                       r = 0, mortality_regime = NULL) { # mortality regime must be specified if age_structured = TRUE
  if(age_structured == TRUE){
    if(is.null(mortality_regime)){
      print("You need to specify a mortality regime in order to generate an age-structured population. Check your function arguments. Does mortality_regime = NULL?")
    }
    pop0 <- create_pop_stable_age(pop0_size = pop0_size,
                                  mortality_regime = mortality_regime,
                                  r       = r, # pop. growth rate
                                  max_age = 100,
                                  year_born = if_else(age == 0, 0, NA_real_)
                                  )
  }
  
  if(age_structured == FALSE){
    pop0 <- data.frame(agent_id = 1:pop0_size,
                       age = 0
    )
  }
  if(pop_config$model_lesions){
    pop0 <- pop0 %>%
      mutate(lesion = if_else(pop0$age %in% pop_config$lesion_formation_window[1]:pop_config$lesion_formation_window[2], 0, NA_real_)) %>% 
      relocate(lesion, .after = age) # change position of lesion column so it sits to the right of 'age'
  }
  if(pop_config$model_frailty){
    if(is.null(c(pop_config$gammafrailty_shape, pop_config$gammafrailty_scale))){
      print("model_frailty = TRUE, but you have not specified shape and scale arguments for the frailty distribution.")
    }
    pop0$frailty <- rgamma(pop0_size, shape = pop_config$gammafrailty_shape, scale = pop_config$gammafrailty_scale)
    pop0$acquired_frailty = NA_real_
  }
  if (!is.null(pop_config$annual_exposure)) {
    pop0$n_stress_events <- 0L
  }
  return(pop0)
}




