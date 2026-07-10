#' Resolve which agents die this step, assign cause, and split the population
#'
#' Given illness and trauma hazards for each agent, converts them to a
#' death probability, draws the death dice, optionally overrides everyone
#' to die (\code{force_death}), assigns cause of death among those who died
#' (proportional to relative cause-specific hazard), and returns the
#' population split into survivors and decedents with step-level columns
#' cleaned up appropriately.
#'
#' This function deliberately does NOT know how illness_hazard/trauma_hazard
#' were computed -- it takes them as given and answers exactly one
#' question: "who dies, from what, and what does the population look like
#' afterward?" Keeping this isolated makes the death-dice and
#' cause-assignment logic testable independent of the hazard-assembly
#' machinery upstream.
#'
#' ## Cause assignment
#' Standard competing-risks cause assignment: among those who died this
#' step, cause is assigned with probability proportional to each agent's
#' relative cause-specific hazard at the moment of death, via a second,
#' independent uniform draw compared against the illness/trauma hazard
#' ratio. This is NOT a second all-or-nothing dice roll per cause, which
#' would double-count risk and inflate total mortality above the death
#' probability implied by illness_hazard + trauma_hazard.
#'
#' ## force_death
#' When \code{force_death = TRUE}, every agent's death probability is set
#' to 1, but illness_hazard/trauma_hazard are still used, unmodified, for
#' cause assignment -- so agents forced to die (e.g. in a population-
#' truncation step) still get a cause of death reflecting the real
#' age-specific illness/trauma split, not an arbitrary default.
#'
#' ## cause_of_death column behavior
#' If \code{split_active = FALSE} (no cause-of-death splitting in effect
#' this step -- see \code{compute_cause_specific_hazards()}), no
#' \code{cause_of_death} column is added to decedents at all, preserving
#' exact backward compatibility with pre-competing-hazards decedent data
#' frames. This matters for \code{do.call(rbind, decedents)} in
#' \code{Simulate_Cemetery}: every batch of decedents across the whole run
#' must have identical columns, so \code{split_active} must be consistent
#' across steps (it will be, since it depends only on \code{mu0_lookup},
#' which is fixed for the whole simulation run).
#'
#' @param pop Population data frame.
#' @param illness_hazard Numeric vector, from
#'   \code{compute_cause_specific_hazards()}.
#' @param trauma_hazard Numeric vector, from
#'   \code{compute_cause_specific_hazards()}.
#' @param split_active Logical scalar, from
#'   \code{compute_cause_specific_hazards()}.
#' @param current_time Integer, current time step (recorded as
#'   \code{year_died} on decedents).
#' @param force_death Logical. If TRUE, every agent dies this step
#'   regardless of hazard (cause assignment still uses the real hazards).
#'
#' @return A list with:
#'   \item{pop}{Survivors, with step-level columns stripped as appropriate.}
#'   \item{decedents}{Decedents this step, with \code{year_died} (and
#'     \code{cause_of_death}, if \code{split_active}) added.}
#' @keywords internal
resolve_deaths <- function(pop,
                           illness_hazard,
                           trauma_hazard,
                           split_active,
                           current_time,
                           force_death              = FALSE) {
  # Combined probability of dying at this time point, from any cause
  threshold <- 1 - exp(-(illness_hazard + trauma_hazard))

  if (force_death) {
    threshold <- rep(1, nrow(pop))
  }

  death_dice <- runif(n = nrow(pop), min = 0, max = 1)
  died <- death_dice < threshold

  # --- Assign cause of death, conditional on death -----------------------
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

  # --- Split into decedents and survivors, clean up step-level columns --
  decedents_this_step <- pop[died, ]
  if (nrow(decedents_this_step) > 0) {
    decedents_this_step$year_died <- current_time
    if (split_active) {
      decedents_this_step$cause_of_death <- cause_of_death[died]
    }
  }

  survivors <- pop[!died, ]

  list(
    pop       = survivors,
    decedents = decedents_this_step
  )
}
