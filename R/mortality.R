#' @title Apply Mortality to Agent Population
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

