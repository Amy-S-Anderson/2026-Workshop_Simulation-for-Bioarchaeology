
#' Sample which agents are exposed to lesion-causing conditions this time step
#'
#' Runs before form_lesions() and apply_mortality() each time step.
#' Writes exposed_this_step (logical) and increments n_stress_events for
#' exposed agents.
#'
#' @param pop Population data frame
#' @param annual_exposure Numeric proportion of agents exposed each time step. Deterministic.
#' @return Updated pop data frame
#' @keywords internal
#' @export
sample_exposure <- function(pop, annual_exposure) {
  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  
  if(annual_exposure < 0 | annual_exposure > 1){
    stop("annual_exposure is the proportion of the population exposed each time step, and must be between 0 and 1")
  }
  
  
  # ---------------------------------------------------------------------------
  # All agents are still unexposed
  # ---------------------------------------------------------------------------
  
  # No one has been exposed yet. After the first time step, this will overwrite exposed_this_step values from the previous time step, so that the output of sample_exposure returns true values for exposure status in the current time step. 
  pop$exposed_this_step <- rep(FALSE, size = nrow(pop))
  
  # ---------------------------------------------------------------------------
  # Identify exposed agents
  # ---------------------------------------------------------------------------
  
  stress <- vector(mode = "logical", length = nrow(pop))
  n_exposed <- round(nrow(pop) * annual_exposure)
  stress[1:n_exposed] <- TRUE
  pop$exposed_this_step <- sample(stress, size = nrow(pop), replace = FALSE)
  
  if ("n_stress_events" %in% names(pop)) {
    pop$n_stress_events <- pop$n_stress_events + as.integer(pop$exposed_this_step)
  }
  
  pop
}
