#' Compute the per-agent hazard multiplier from risk factors and exposure
#'
#' Combines all step-level, non-age-based modifiers to mortality hazard into
#' a single multiplier per agent: general risk factors (e.g. frailty,
#' acquired_frailty, lesion presence) and, if applicable, a transient hazard
#' from exposure to lesion-causing conditions this step.
#'
#' This function deliberately does NOT know about ages, cause-specific
#' hazards, or death. It answers exactly one question: "what number do I
#' multiply this agent's baseline hazard by, this step, because of their
#' traits and exposures?" Keeping that question isolated here makes it
#' possible to unit-test the risk-factor and exposure logic without any of
#' the mortality machinery (death dice, cause assignment, etc.) running.
#'
#' ## Exposure-hazard logic
#' Two independent switches govern whether/how exposure to lesion-causing
#' conditions affects hazard:
#'   - \code{exposure_causes_hazard}: does exposure add mortality hazard at
#'     all, independent of lesion status itself?
#'   - \code{hazard_is_transient}: if so, is that hazard transient (applies
#'     only during the step of exposure, computed here) or permanent
#'     (becomes ongoing acquired frailty, computed elsewhere)?
#'
#' \code{lesion_requires_survival} does NOT gate whether exposure causes
#' hazard -- it only determines, in \code{Simulate_Cemetery}'s time loop,
#' whether lesion formation happens before or after the mortality check
#' each step. This function's exposure-hazard branch applies identically
#' regardless of \code{lesion_requires_survival}, provided
#' \code{exposed_this_step} exists on \code{pop} (which it will, in either
#' ordering, since \code{sample_exposure()} runs before either ordering
#' begins).
#'
#' When \code{hazard_is_transient = FALSE} (permanent exposure hazard),
#' this function does nothing further: permanent hazard is expected to
#' already be reflected in a persistent column (e.g. \code{acquired_frailty}
#' or \code{lesion}) that risk_factors picks up, because \code{form_lesions()}
#' is responsible for converting exposure into that ongoing state. This
#' function only ever computes the *transient* piece directly.
#'
#' @param pop Population data frame. Must contain any columns named in
#'   \code{risk_factors}, and \code{exposed_this_step} if
#'   \code{exposure_causes_hazard = TRUE}.
#' @param risk_factors Named list of risk factor columns. For continuous
#'   columns (e.g. frailty, acquired_frailty) whose values ARE the
#'   proportional hazard, set the value to NULL. For binary (0/1) columns
#'   (e.g. lesion), supply the proportional hazard as a numeric scalar
#'   (e.g. list(lesion = 2.1)).
#' @param exposure_hazard_multiplier Numeric. The multiplier applied to
#'   exposed agents when the exposure hazard is transient.
#'
#' @return Numeric vector, length \code{nrow(pop)}, giving the combined
#'   hazard multiplier for each agent.
#' @keywords internal
compute_hazard_multiplier <- function(pop,
                                      risk_factors               = list(),
                                      exposure_hazard_multiplier = 1) {
  
  hazard_multiplier <- rep(1, nrow(pop))
  
  # --- General risk factors (frailty, acquired_frailty, etc.) ---
  # NULL entries: column value IS the multiplier (e.g. continuous frailty).
  # Numeric scalar entries: binary column; if 1, multiply by this scalar.
  for (factor_name in names(risk_factors)) {
    if (!factor_name %in% names(pop)) next
    col <- pop[[factor_name]]
    col[is.na(col)] <- 1
    factor_spec <- risk_factors[[factor_name]]

    if (is.null(factor_spec)) {
      hazard_multiplier <- hazard_multiplier * col # loops through named risk factors sequentially, with the previous hazard multiplier output now as an input value.
    } else if (is.numeric(factor_spec) && length(factor_spec) == 1) {
      hazard_multiplier <- hazard_multiplier * ifelse(col == 1, factor_spec, 1)
    } else {
      warning(sprintf("risk_factors[['%s']] must be NULL or a single numeric. Skipping.",
                      factor_name))
    }
  }
  

  hazard_multiplier
}