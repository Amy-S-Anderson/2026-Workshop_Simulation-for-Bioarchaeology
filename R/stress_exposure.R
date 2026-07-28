
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
sample_exposure <- function(pop, annual_exposure, 
                            lesion_formation_rate, 
                            lesion_formation_window) {
  
  # ---------------------------------------------------------------------
  # Input validation -- guard with is.null() first, since exactly one of
  # annual_exposure / lesion_formation_rate is expected to be NULL by
  # design (see form_lesions()'s mutual-exclusivity check). Comparing
  # NULL < 0 returns logical(0), and if(logical(0)) errors -- so these
  # checks must not run at all when the argument is NULL.
  # ---------------------------------------------------------------------
  if (!is.null(annual_exposure) && (annual_exposure < 0 || annual_exposure > 1)) {
    stop("annual_exposure is the proportion of the population exposed each time step, and must be between 0 and 1")
  }
  if (!is.null(lesion_formation_rate) && (lesion_formation_rate < 0 || lesion_formation_rate > 1)) {
    stop("lesion_formation_rate is the proportion of agents in the lesion formation age window who are exposed to lesion-causing conditions each time step, and must be between 0 and 1")
  }
  
  # ---------------------------------------------------------------------
  # All agents start unexposed this step
  # ---------------------------------------------------------------------
  pop$exposed_this_step <- FALSE  # recycles to nrow(pop) automatically
  
  # ---------------------------------------------------------------------
  # Identify exposed agents
  # ---------------------------------------------------------------------
  in_window <- pop$age >= lesion_formation_window[1] &
    pop$age <= lesion_formation_window[2]
  
  if (!is.null(annual_exposure)) {  # annual exposure framework
    n_exposed <- round(nrow(pop) * annual_exposure)
    stress <- rep(FALSE, nrow(pop))
    if (n_exposed > 0) stress[seq_len(n_exposed)] <- TRUE
    pop$exposed_this_step <- sample(stress, size = nrow(pop), replace = FALSE)
    
  } else {  # lesion formation framework
    n_in_window <- sum(in_window)
    if (n_in_window > 0) {
      n_exposed <- round(n_in_window * lesion_formation_rate)
      stress <- rep(FALSE, n_in_window)
      if (n_exposed > 0) stress[seq_len(n_exposed)] <- TRUE
      pop$exposed_this_step[in_window] <- sample(stress, size = n_in_window, replace = FALSE)
    }
  }
  
  if ("n_stress_events" %in% names(pop)) {
    pop$n_stress_events <- pop$n_stress_events + as.integer(pop$exposed_this_step)
  }
  
  pop
}

