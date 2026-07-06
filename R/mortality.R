#' Helper functions for calculating mortality risk

#... I think these might be unnecessary now that defrail_siler is in the mix. 


#'@title Siler hazard (helper)
#'@description Compute Siler hazard for a given age. This helper function is referenced in apply_mortality(). 
#'
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


#' @title Siler survivorship (helper)
#' @description Compute discrete Siler survivorship l(x). Returns the probability of surviving from birth to each age x, under the Siler hazard model. l(0) is defined as 1 (everyone is alive at birth).
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




#' @title Apply Mortality to Agent Population
#'
#' @description Applies a single time step of mortality to a population of
#'   agents, based on age-dependent mortality hazards (defined by Siler
#'   function parameters) and, optionally, other sources of mortality risk
#'   (e.g., individual frailty, skeletal lesion presence, etc.). Optionally
#'   splits mortality into competing "illness" and "trauma" causes, so that
#'   risk factors (e.g., lesions) can be specified to affect illness-related
#'   mortality without affecting trauma-related mortality.
#'
#' @param pop Population data frame
#' @param current_time Integer, the current time step in the model
#' @param mu0_lookup data.frame, created by defrail_siler at the top of
#'   Simulate_Cemetery, optionally merged with a p_trauma column via
#'   build_p_trauma_lookup(). Defaults to NULL, in which case there is no
#'   heterogeneity in individual frailty at birth. If mu0_lookup has no
#'   p_trauma column, or p_trauma is 0 at every age, all mortality is
#'   attributed to illness and no cause_of_death column is added to
#'   decedents (fully backward compatible behavior).
#' @param schedule Named vector of Siler parameter values. The mortality
#'   regime for baseline hazards, before being modified by individual
#'   traits/states beyond age.
#' @param mortality_risk_type One of "proportional", "time_decreasing",
#'   "time_increasing"
#' @param force_death Logical. Defaults to FALSE. If TRUE, kill everyone at this time step. 
#' @param risk_factors Named list of risk factor columns. For continuous
#'   columns (e.g. frailty, acquired_frailty) whose values ARE the
#'   proportional hazard, set the value to NULL. For binary (0/1) columns
#'   (e.g. lesion), supply the proportional hazard as a numeric scalar
#'   (e.g. list(lesion = 2.1)). NOTE: when cause-of-death splitting is
#'   active (i.e. mu0_lookup$p_trauma has nonzero values), all risk_factors
#'   are currently applied to illness-related mortality only; trauma-related
#'   mortality reflects the age-based baseline hazard alone. This matches
#'   the intended use case (e.g. skeletal lesions raising risk of dying from
#'   illness/infection but not from accidents), but means risk_factors
#'   intended to modify trauma risk are not yet supported.
#' @return Updated pop data frame
#' @keywords internal
#' @export

apply_mortality <- function(pop,
                            mu0_lookup,
                            mortality_risk_type         = "proportional",
                            force_death                 = FALSE,
                            risk_factors                = list(),
                            current_time,
                            lesion_requires_survival    = FALSE,
                            exposure_causes_hazard      = FALSE,
                            hazard_is_transient         = FALSE,
                            exposure_hazard_multiplier  = 1,
                            lesion_formation_window     = NULL) {
  
  #
  age_based_risk <- mu0_lookup$mu0[match(pop$age, mu0_lookup$age)]
  
  
  # --- Competing hazards: split age-based risk into illness/trauma ---------
  
  # If mu0_lookup carries no p_trauma column (or it is entirely 0/NA), the
  # split proportion is 0 everywhere, illness gets the full age-based risk,
  # and trauma is 0. This means any lesion-related mortality risk is applied to everyone, and all-cause mortality functions as the single mortality category. 
  # if mu0_lookup does have a p_trauma column, then use those age-specific probabilities of traumatic death (which will be unaffected by lesion-causing conditions/lesion status)
  if ("p_trauma" %in% names(mu0_lookup)) {
    p_trauma_by_age <- mu0_lookup$p_trauma[match(pop$age, mu0_lookup$age)]
    p_trauma_by_age[is.na(p_trauma_by_age)] <- 0
  } else {
    p_trauma_by_age <- rep(0, nrow(pop))
  }
  # Indicate whether cause-of-death splits are being tracked in this model
  split_active <- any(p_trauma_by_age > 0, na.rm = TRUE)
  
  # If they are, then calculate probability of death for: 
  age_based_illness <- age_based_risk * (1 - p_trauma_by_age)
  age_based_trauma  <- age_based_risk * p_trauma_by_age
  
  # fallback hazard multiplier = 1 (no additional hazards)
  hazard_multiplier <- 1
  
  for (factor_name in names(risk_factors)) { # right now, risk_factors are things like 'frailty' and 'lesion'
    if (!factor_name %in% names(pop)) next
    col <- pop[[factor_name]] # frailty value, or lesion-related risk
    col[is.na(col)] <- 1
    factor_spec <- risk_factors[[factor_name]]
    
    # if no risk factors are present
    if (is.null(factor_spec)) {
      hazard_multiplier <- hazard_multiplier * col # still 1
    } else if (is.numeric(factor_spec) && length(factor_spec) == 1) {
      hazard_multiplier <- hazard_multiplier * ifelse(col == 1, factor_spec, 1) 
    } else {
      warning(sprintf("risk_factors[['%s']] must be NULL or a single numeric. Skipping.",
                      factor_name))
    }
  }
  
  # Apply transient hazard if present
  if ("transient_hazard" %in% names(pop)) {
    hazard_multiplier <- hazard_multiplier * pop$transient_hazard
  }
  
  # In lesion_requires_survival cases, apply exposure hazard before mortality check
  # (form_lesions hasn't run yet, so we apply the hazard directly here)
  if (lesion_requires_survival && exposure_causes_hazard && "exposed_this_step" %in% names(pop)) {
    if (hazard_is_transient) {
      hazard_multiplier <- hazard_multiplier * ifelse(pop$exposed_this_step,
                                                      exposure_hazard_multiplier, 1)
    } else {
      # Permanent frailty for lesion_requires_survival cases is handled in
      # form_lesions() after this function returns, for survivors only
    }
  }
  
  # --- Assemble cause-specific hazards ---------------------------------
  # risk_factors, pop$frailty, and hazard_multiplier (lesion/exposure/
  # transient effects) modify illness-related mortality only. Trauma-related
  # mortality reflects the age-based baseline hazard (age_based_trauma)
  # alone, so that e.g. lesion status can raise risk of dying from illness
  # without raising risk of dying from trauma.
  if (mortality_risk_type == "proportional") {
    illness_hazard <- age_based_illness * pop$frailty * hazard_multiplier
    trauma_hazard  <- age_based_trauma
  } else if (mortality_risk_type == "time_decreasing") {
    illness_hazard <- age_based_illness * pop$frailty * hazard_multiplier /
      ((pop$age / 10) + hazard_multiplier)
    trauma_hazard  <- age_based_trauma
  } else if (mortality_risk_type == "time_increasing") {
    illness_hazard <- age_based_illness * pop$frailty *
      ((pop$age / 10) + hazard_multiplier) / hazard_multiplier
    trauma_hazard  <- age_based_trauma
  } else {
    stop(sprintf("Unrecognized mortality_risk_type: '%s'", mortality_risk_type))
  }
  
  threshold <- 1 - exp(-(illness_hazard + trauma_hazard))
  
  # if force_death = TRUE, kill everyone. 
  if (force_death) {
    threshold <- rep(1, nrow(pop))
  }
  
  death_dice <- runif(n = nrow(pop), min = 0, max = 1)
  died <- death_dice < threshold
  
  # --- Assign cause of death, conditional on death ----------------------
  # Standard competing-risks cause assignment: among those who died this
  # step, cause is assigned with probability proportional to each agent's
  # relative cause-specific hazard at the moment of death. This is done
  # with a second, independent uniform draw compared against the illness/
  # trauma hazard ratio -- NOT a second all-or-nothing dice roll per cause,
  # which would double-count risk and inflate total mortality above
  # `threshold`.
  cause_of_death <- NULL
  if (split_active) {
    cause_of_death <- rep(NA_character_, nrow(pop))
    if (any(died)) {
      total_hazard_died <- illness_hazard[died] + trauma_hazard[died]
      # Guard against 0/0 (can only occur if threshold was 0 for a death,
      # which shouldn't happen, but default such cases to illness rather
      # than propagating NaN)
      p_illness_given_death <- ifelse(
        total_hazard_died > 0,
        illness_hazard[died] / total_hazard_died,
        1
      )
      cause_dice <- runif(n = sum(died), min = 0, max = 1)
      cause_of_death[died] <- ifelse(cause_dice < p_illness_given_death,
                                     "illness", "trauma")
    }
  }
  
  # Strip step-level columns from decedents unconditionally
  decedents_this_step <- pop[died, ]
  if ("transient_hazard"   %in% names(decedents_this_step)) decedents_this_step$transient_hazard   <- NULL
  if ("exposed_this_step"  %in% names(decedents_this_step)) decedents_this_step$exposed_this_step  <- NULL
  if (nrow(decedents_this_step) > 0) {
    decedents_this_step$year_died <- current_time
    if (split_active) {
      decedents_this_step$cause_of_death <- cause_of_death[died]
    }
  }
  
  # Strip step-level columns from survivors, EXCEPT keep exposed_this_step
  # when lesion_requires_survival = TRUE so form_lesions() can still read it
  survivors <- pop[!died, ]
  if ("transient_hazard" %in% names(survivors)) survivors$transient_hazard <- NULL
  if (!lesion_requires_survival && "exposed_this_step" %in% names(survivors)) {
    survivors$exposed_this_step <- NULL
  }
  
  list(
    pop       = survivors,
    decedents = decedents_this_step
  )
}


#' 
#' 
#' 
#' 
#' 
#' 
#' @title Apply mortality using usher3
#'
#' @description Applies a single time step of mortality to a population of
#'   agents using the Usher3 illness-death model. Agents with lesions
#'   experience elevated mortality risk compared to agents without lesions.
#'
#' @details This function implements a discrete-time approximation of the
#'   Usher3 illness-death model for mortality. Within each time step of size
#'   \code{dx}, each living agent faces a probability of death equal to
#'   \code{hazard * dx}, where the hazard depends on the agent's lesion status:
#'
#'   \itemize{
#'     \item Agents without lesions: \code{p_die = hsiler(age, b_siler) * dx}
#'     \item Agents with lesions: \code{p_die = k2 * hsiler(age, b_siler) * dx}
#'   }
#'
#'   This approximation is valid when \code{hazard * dx << 1}. For accuracy,
#'   use small values of \code{dx} (e.g., 0.01 to 1).
#'
#'   Agents who die retain their age at death. Surviving agents have their
#'   age incremented by \code{dx}.
#'
#'   Note: This module handles mortality only. Lesion acquisition and age
#'   filtration (archaeological sampling bias) are handled by separate modules.
#'
#' @param state A data.frame representing the current population state. Must
#'   contain the following columns:
#'   \describe{
#'     \item{agent_id}{Integer. Unique identifier for each agent.}
#'     \item{age}{Numeric. Current age of each agent.}
#'     \item{lesion}{Logical. TRUE if the agent has a lesion, FALSE otherwise.}
#'     \item{dead}{Logical. TRUE if the agent is dead, FALSE if alive.}
#'     \item{in_sample}{Logical. TRUE if the agent is in the sample.}
#'   }
#'
#' @param mortality_model Character. The mortality model to use. Currently
#'   only \code{"usher3"} is supported. Any other value will raise an error.
#'
#' @param mortality_param Numeric vector containing the mortality parameters.
#'   For the \code{"usher3"} model, this must be a length-6 vector:
#'   \describe{
#'     \item{\code{[1]} k2}{The mortality multiplier for agents with lesions. A value
#'       of 1.0 means no excess mortality; values > 1 indicate elevated risk.
#'       Must be non-negative.}
#'     \item{\code{[2:6]} b_siler}{The five Siler hazard parameters using the demohaz
#'       parameterization: \code{c(b1, b2, b3, b4, b5)}. See
#'       \code{\link[demohaz]{hsiler}} for details.}
#'   }
#'
#' @param dx Numeric. The time step size in years. Determines both the
#'   mortality probability (hazard * dx) and the age increment for survivors.
#'   Smaller values give more accurate results but require more iterations.
#'
#' @return A data.frame with the same structure as \code{state}, with updated
#'   values:
#'   \itemize{
#'     \item \code{dead}: Updated to TRUE for agents who died this time step.
#'     \item \code{age}: Incremented by \code{dx} for surviving agents;
#'       unchanged for agents who died (preserving age at death).
#'     \item \code{in_sample}: Unchanged (passed through from input).
#'   }
#'
#' @examples
#' # Create a small test population
# state <- data.frame(
#   agent_id = 1:10,
#   age = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80),
#   lesion = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE),
#   dead = rep(FALSE, 10),
#   in_sample = rep(TRUE, 10)
# )
# 
# # Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
# a_siler <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
# b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
# mortality_param <- c(1.2, b_siler)  # k2 = 1.2
# 
# # Apply one year of mortality
# result <- apply_mortality_usher3(
#   state = state,
#   mortality_model = "usher3",
#   mortality_param = mortality_param,
#   dx = 1
# )
#'
#' @seealso \code{\link[demohaz]{hsiler}} for the Siler hazard function.
#'
#' @export
apply_mortality_usher3 <- function(state,
                            mortality_model,
                            mortality_param,
                            dx) {


  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (is.null(mortality_model)) {
    stop("mortality_model must be specified")
  }

  # ---------------------------------------------------------------------------
  # Extract parameters (model-specific)
  # ---------------------------------------------------------------------------

  if (mortality_model == "usher3") {
    if (length(mortality_param) != 6) {
      stop("For usher3, mortality_param must be length 6: c(k2, b1, b2, b3, b4, b5)")
    }
    k2 <- mortality_param[1]
    b_siler <- mortality_param[2:6]
  } else {
    stop("mortality_model must be 'usher3'. Other models are not yet supported.")
  }

  # ---------------------------------------------------------------------------
  # Copy state to avoid modifying input
  # ---------------------------------------------------------------------------

  result <- state

  # ---------------------------------------------------------------------------
  # Identify living agents
  # ---------------------------------------------------------------------------

  alive_idx <- which(!state$dead)

  # If no one is alive, return early
  if (length(alive_idx) == 0) {
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Calculate death probabilities for living agents
  # ---------------------------------------------------------------------------

  # Current ages and lesion status of living agents
  x <- state$age[alive_idx]
  has_lesion <- state$lesion[alive_idx]

  # Baseline Siler mortality hazard at each agent's current age

  h <- demohaz::hsiler(x, b_siler)

  # Death probability depends on lesion status:
  #   - Without lesion: p_die = h * dx

  #   - With lesion:    p_die = k2 * h * dx
  p_die <- ifelse(has_lesion, k2 * h * dx, h * dx)

  # ---------------------------------------------------------------------------
  # Determine who dies this time step
  # ---------------------------------------------------------------------------

  # One uniform draw per living agent

  u <- runif(length(alive_idx))

  # Agent dies if their random draw is below their death probability

  dies <- u < p_die

  # Indices in the original state data.frame

  dead_idx <- alive_idx[dies]
  survivor_idx <- alive_idx[!dies]

  # ---------------------------------------------------------------------------
  # Update state
  # ---------------------------------------------------------------------------

  # Mark deaths
  result$dead[dead_idx] <- TRUE

  # Increment age for survivors only

  # Dead agents retain their age at death (no increment)
  result$age[survivor_idx] <- result$age[survivor_idx] + dx

  # ---------------------------------------------------------------------------
  # Return updated state
  # ---------------------------------------------------------------------------

  return(result)
}

