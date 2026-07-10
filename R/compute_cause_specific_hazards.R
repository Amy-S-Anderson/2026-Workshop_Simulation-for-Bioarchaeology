#' Compute cause-specific (illness vs. trauma) hazards for each agent
#'
#' Given each agent's age-based baseline risk, frailty, and hazard
#' multiplier (from \code{compute_hazard_multiplier()}), splits mortality
#' hazard into illness and trauma components using the age-specific
#' \code{p_trauma} proportions in \code{mu0_lookup}, and applies the chosen
#' \code{mortality_risk_type} functional form.
#'
#' This function deliberately does NOT know about risk_factors, exposure,
#' or death dice -- it takes the hazard_multiplier as a given input and
#' answers exactly one question: "what are this agent's illness and trauma
#' hazards this step?" Keeping this isolated makes the age/cause-splitting
#' math (and the mortality_risk_type branches) testable independent of
#' how hazard_multiplier was assembled.
#'
#' ## Cause splitting
#' If \code{mu0_lookup} has no \code{p_trauma} column, or \code{p_trauma}
#' is 0 at every age present in \code{pop}, all hazard is attributed to
#' illness and \code{trauma_hazard} is 0 everywhere -- this reproduces
#' pre-competing-hazards behavior exactly (see \code{split_active} in the
#' return value).
#'
#' ## Frailty/hazard_multiplier only affects illness
#' By design (per the project's modeling goals), \code{pop$frailty} and
#' \code{hazard_multiplier} (which carries lesion/exposure/transient
#' effects) modify illness-related hazard only. Trauma-related hazard is
#' the age-based baseline alone, so that e.g. lesion status can raise risk
#' of dying from illness without raising risk of dying from trauma. This is
#' NOT a general property of competing-hazards models -- it is a deliberate
#' simplification for this project's current research question. If future
#' work needs risk factors that modify trauma hazard too (e.g. sex effects
#' on accident risk), this function will need a parallel
#' \code{trauma_risk_factors} pathway; it is not there yet.
#'
#' @param pop Population data frame. Must contain \code{age} and
#'   \code{frailty} columns.
#' @param mu0_lookup Data frame with \code{age} and \code{mu0} columns
#'   (from \code{defrail_siler()}), optionally with a \code{p_trauma}
#'   column (from \code{build_p_trauma_lookup()}).
#' @param hazard_multiplier Numeric vector, length \code{nrow(pop)}, from
#'   \code{compute_hazard_multiplier()}.
#' @param mortality_risk_type One of "proportional", "time_decreasing",
#'   "time_increasing".
#'
#' @return A list with:
#'   \item{illness_hazard}{Numeric vector, illness-attributable hazard per agent.}
#'   \item{trauma_hazard}{Numeric vector, trauma-attributable hazard per agent.}
#'   \item{split_active}{Logical scalar. TRUE if any agent has nonzero
#'     p_trauma, i.e. cause-of-death splitting is meaningfully in effect
#'     this step.}
#' @keywords internal
compute_cause_specific_hazards <- function(pop,
                                           mu0_lookup,
                                           hazard_multiplier,
                                           mortality_risk_type = "proportional") {

  age_based_risk <- mu0_lookup$mu0[match(pop$age, mu0_lookup$age)]

  # --- Split age-based risk into illness/trauma proportions -------------
  if ("p_trauma" %in% names(mu0_lookup)) {
    p_trauma_by_age <- mu0_lookup$p_trauma[match(pop$age, mu0_lookup$age)]
    p_trauma_by_age[is.na(p_trauma_by_age)] <- 0
  } else {
    p_trauma_by_age <- rep(0, nrow(pop))
  }
  split_active <- "p_trauma" %in% names(mu0_lookup) && any(mu0_lookup$p_trauma > 0, na.rm = TRUE)
  
  age_based_illness <- age_based_risk * (1 - p_trauma_by_age)
  age_based_trauma  <- age_based_risk * p_trauma_by_age

  # --- Apply mortality_risk_type functional form -------------------------
  # frailty and hazard_multiplier modify illness_hazard only; trauma_hazard
  # is the age-based baseline alone (see function-level doc).
  if (mortality_risk_type == "proportional") {
    illness_hazard <- age_based_illness * hazard_multiplier # hazard_multiplier contains frailty
    trauma_hazard  <- age_based_trauma
  } else if (mortality_risk_type == "time_decreasing") {
    illness_hazard <- age_based_illness  * hazard_multiplier /
      ((pop$age / 10) + hazard_multiplier)
    trauma_hazard  <- age_based_trauma
  } else if (mortality_risk_type == "time_increasing") {
    illness_hazard <- age_based_illness * 
      ((pop$age / 10) + hazard_multiplier) / hazard_multiplier
    trauma_hazard  <- age_based_trauma
  } else {
    stop(sprintf("Unrecognized mortality_risk_type: '%s'", mortality_risk_type))
  }

  list(
    illness_hazard = illness_hazard,
    trauma_hazard  = trauma_hazard,
    split_active   = split_active
  )
}
