


#'@title Update acquired frailty
#'@description Compute new acquired frailty as a function of exposure to stress.
#'
#' @param ages Integer vector, if population is age-structured
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric hazard value
#' @keywords internal
#' 

acquire_frailty <- function(pop, exposure_hazard_multiplier){
  newly_acquired_frailty <- ifelse(pop$exposed_this_step, exposure_hazard_multiplier, NA_real_)
  old_acquired_frailty <- pop$acquired_frailty
  updated_frailty <- ifelse()
  
  pop$acquired_frailty <- updated_frailty
}


#' Update acquired_frailty for one time step of exposure
#'
#' Computes this step's newly acquired frailty (if the agent was exposed
#' this step) and adds it to the agent's existing accumulated frailty.
#'
#' NA represents "has not yet acquired any frailty" -- kept distinct from
#' a literal 0 so that acquired_frailty's "no effect yet" sentinel never
#' collides with 0 as a legitimately-coded value in some OTHER risk_factors
#' column (e.g. lesion, where 0 means "confirmed no lesion", not "unknown").
#' See compute_hazard_multiplier() for the generic NA -> 1 handling this
#' choice enables safely across every risk factor column.
#'
#' NA-aware addition: if both old and new values are NA (never exposed,
#' this step or previously), the result is NA. Otherwise NA is treated as
#' 0 and the two values are summed.
#'
#' @param pop Population data frame. Must contain acquired_frailty and
#'   exposed_this_step columns.
#' @param exposure_hazard_multiplier Numeric. Frailty increment added for
#'   a single exposure event this step.
#' @return pop, with acquired_frailty updated. Caller MUST use this
#'   returned pop going forward (same rule as compute_hazard_multiplier's
#'   pop -- see that function's docs for what happens if you don't).
#' @keywords internal
acquire_frailty <- function(pop, exposure_hazard_multiplier) {
  
  newly_acquired_frailty <- ifelse(pop$exposed_this_step,
                                   exposure_hazard_multiplier,
                                   NA_real_)
  old_acquired_frailty <- pop$acquired_frailty
  
  both_na <- is.na(old_acquired_frailty) & is.na(newly_acquired_frailty)
  summed  <- dplyr::coalesce(old_acquired_frailty, 0) + # coalesce replaces NA values with a specified value, or with non-missing values from the second element/vector in the function call. 
    dplyr::coalesce(newly_acquired_frailty, 0) # 
  
  pop$acquired_frailty <- ifelse(both_na, NA_real_, summed)
  pop
}