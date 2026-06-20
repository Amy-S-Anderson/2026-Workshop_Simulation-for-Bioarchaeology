#' @title Apply Deposition Filtering to Agent Population
#'
#' @description Applies deposition filtering to a population of dead agents.
#'   Determines which deceased individuals are deposited into the
#'   archaeological record based on age cutoff criteria.
#'
#' @details This function implements the first stage of archaeological
#'   filtration. Dead agents whose age at death meets the deposition criteria
#'   are marked as deposited (\code{was_deposited = TRUE}).
#'
#'   Currently, only the \code{"cutoff"} model is supported, which deposits
#'   all dead agents at or above a minimum age. This models scenarios where
#'   very young individuals (e.g., neonates) are less likely to enter the
#'   archaeological record due to burial practices or taphonomic factors.
#'
#'   The deposition decision is deterministic: all dead agents at or above
#'   the cutoff age are deposited; all below are not.
#'
#'   Note: This module handles deposition only and is applied to the decedents data frame. Living agents are not
#'   processed. Preservation filtering (whether deposited remains are
#'   recoverable) is handled by a separate module.
#'
#' @param decedents A data.frame representing the full dead population.
#'   will contain the following columns, and possibly some other agent states:
#'   \describe{
#'     \item{agent_id}{Integer. Unique identifier for each agent.}
#'     \item{age}{Numeric. Age at death for dead agents.}
#'     \item{was_deposited}{Logical. Will be updated by this function.}
#'     \item{in_sample}{Logical. TRUE if the agent is in the sample.}
#'   }
#'
#' @param deposition_model Character. The deposition model to use. Currently
#'   only \code{"cutoff"} is supported.
#'
#' @param deposition_param Numeric vector containing the deposition parameters.
#'   For the \code{"cutoff"} model, this must be length 1:
#'   \describe{
#'     \item{cutoff_age}{The minimum age for deposition. Dead agents with
#'       age >= cutoff_age are deposited; those below are not.}
#'   }
#'
#' @param dx Numeric. The time step size in years. Not used in the current
#'   implementation but included for interface consistency with other modules.
#'
#' @return A data.frame with the same structure as \code{state}, with updated
#'   values:
#'   \itemize{
#'     \item \code{was_deposited}: Set to TRUE for dead agents at or above
#'       the cutoff age; FALSE otherwise.
#'     \item All other columns: Unchanged.
#'   }
#'
#' @examples
#' # Create a test population of deceased individuals
#' state <- data.frame(
#'   agent_id = 1:10,
#'   age = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7),
#'   lesion = rep(FALSE, 10),
#'   was_deposited = rep(FALSE, 10),
#'   in_sample = rep(TRUE, 10)
#' )
#'
#' # Apply deposition with cutoff at age 2.5 (exclude infants/toddlers)
#' result <- apply_deposition(
#'   state = state,
#'   deposition_model = "cutoff",
#'   deposition_param = c(2.5),
#'   dx = 1
#' )
#'
#' @seealso \code{\link{apply_preservation}} for the preservation filtering
#'   module that follows deposition.
#'
#' @export
apply_deposition <- function(decedents,
                             deposition_model,
                             deposition_param,
                             dx) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (is.null(deposition_model)) {
    stop("deposition_model must be specified")
  }

  # ---------------------------------------------------------------------------
  # Extract parameters (model-specific)
  # ---------------------------------------------------------------------------

  if (deposition_model == "cutoff") {
    if (length(deposition_param) != 1) {
      stop("For cutoff model, deposition_param must be length 1: c(cutoff_age)")
    }
    cutoff_age <- deposition_param[1]

  } else {
    stop("deposition_model must be 'cutoff'. Other models are not yet supported.")
  }

  # ---------------------------------------------------------------------------
  # Copy state to avoid modifying input
  # ---------------------------------------------------------------------------

  result <- decedents

  # ---------------------------------------------------------------------------
  # Identify dead agents eligible for deposition
  # ---------------------------------------------------------------------------

  #.... I'm not sure this if statement needs to be included, now that the function runs on a data frame of only deceased agents, rather than a data frame of all agents with a column to indicate agent state (alive/dead)
  # Only dead agents can be deposited
  # If no one is dead, return early
  # if (nrow(decedents) < 1) {
  #   return(result)
  # }

  # ---------------------------------------------------------------------------
  # Apply deposition logic (model-specific)
  # ---------------------------------------------------------------------------

  if (deposition_model == "cutoff") {
    # Deterministic: deposit if age >= cutoff
    # Mark deposited agents
    result$was_deposited <- decedents$age >= cutoff_age
  }

  # ---------------------------------------------------------------------------
  # Return updated state
  # ---------------------------------------------------------------------------

  return(result)
}

