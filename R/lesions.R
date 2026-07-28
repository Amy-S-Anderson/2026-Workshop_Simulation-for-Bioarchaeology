


#' @title Apply lesion formation before or after mortality
#'
#' @description a single vectorized draw to update each agent's lesion state. Agents who do not have a lesion may acquire one during this action, if they meet the conditional requirements defined by the function arguments. Lesions can either form probabilistically, according to a lesion formation rate, or they can form deterministically, according to the proportion of agents exposed to lesion-causing conditions (annual_exposure). Lesion formation rate will be a better predictor of the number of lesions that form, since annual exposure is spread evenly across agents of all ages, and not all ages are capable of forming lesions (ages of lesion formation are determined by the values of the lesion formation window.)
#'
#'#' Note: These aren't quite parallel options. if(lesion_formation_rate) is used, then acquired frailty values aren't updated. Right now frailty is only integrated with deterministic lesion formation, not probabilistic lesion formation. 
#' @param pop Population data frame
#' @param lesion_formation_window numeric vector of length = 2. c(Age at which lesions can start forming, Age at which lesions stop forming)
#' @param lesion_formation_rate A probability
#' @param annual_exposure A proportion. If lesion_formation_rate has a value then annual_exposure should be set to NULL, and vice versa. 
#' @return Updated pop data frame
#' @keywords internal
#' 
#' 
### Note: Right now, form_lesions doesn't work with lesion_formation_rate, only with annual_exposure. Need to address this. 
#' @export


form_lesions <- function(pop,
                         lesion_formation_window,
                         lesion_formation_rate = NULL,
                         annual_exposure = NULL) {
  
  # -----------------------------------------------------------------------------
  # Input validation
  # -----------------------------------------------------------------------------
  
  # Exactly one of lesion_formation_rate / annual_exposure must be supplied
  n_specified <- sum(!is.null(lesion_formation_rate), !is.null(annual_exposure))
  if (n_specified != 1) {
    stop("Exactly one of lesion_formation_rate or annual_exposure must be specified.")
  }
  
  # Range-check whichever one was supplied
  if (!is.null(lesion_formation_rate) &&
      (lesion_formation_rate < 0 || lesion_formation_rate > 1)) {
    stop("lesion_formation_rate must be between 0 and 1")
  }
  if (!is.null(annual_exposure) &&
      (annual_exposure < 0 || annual_exposure > 1)) {
    stop("annual_exposure must be between 0 and 1")
  }
  
  # lesion_formation_window bounds
  if (lesion_formation_window[1] < 0) {
    stop("formation_window_opens must be non-negative")
  }
  if (lesion_formation_window[2] < lesion_formation_window[1]) {
    stop("formation_window_closes must be >= formation_window_opens")
  }
  
  # -----------------------------------------------------------------------------
  # Identify agents at risk of forming lesions
  # -----------------------------------------------------------------------------
  
  # Identify agents whose ages fall within the range of the lesion formation window. 
  in_window <- pop$age >= lesion_formation_window[1] &
    pop$age <= lesion_formation_window[2]
  
  # exposed_this_step was calculated in sample_exposure()
    exposed <- pop$exposed_this_step
    
    # Form lesions
    pop$lesion <- pmax(pop$lesion, as.integer(in_window & exposed), na.rm = TRUE)
     
     # return the updated population. 
   pop
}









#' @title Apply Lesion Acquisition to Agent Population
#'
#' @description Applies a single time step of lesion acquisition to a
#'   population of agents. Agents without lesions may acquire one based on
#'   a constant transition rate, optionally restricted to an age window.
#'
#' @details This function implements a discrete-time approximation for lesion
#'   acquisition. Within each time step of size \code{dx}, each eligible agent
#'   faces a probability of acquiring a lesion equal to \code{k1 * dx}, where
#'   \code{k1} is the transition rate (hazard coefficient).
#'
#'   An agent is eligible for lesion acquisition if:
#'   \itemize{
#'     \item The agent is alive (\code{dead == FALSE})
#'     \item The agent does not already have a lesion (\code{lesion == FALSE})
#'     \item The agent's age falls within the model's age window
#'   }
#'
#'   This approximation is valid when \code{k1 * dx << 1}. For accuracy,
#'   use small values of \code{dx} (e.g., 0.01 to 1).
#'
#'   Note: This module handles lesion acquisition only. Mortality and age
#'   increments are handled by separate modules.
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
#' @param lesion_model Character. The lesion acquisition model to use. One of:
#'   \describe{
#'     \item{\code{"constant"}}{Constant rate, always possible (no age restriction)}
#'     \item{\code{"constant_to"}}{Constant rate up to a specified age}
#'     \item{\code{"constant_from"}}{Constant rate after a specified age}
#'     \item{\code{"constant_interval"}}{Constant rate within an age interval}
#'   }
#'
#' @param lesion_param Numeric vector containing the lesion parameters.
#'   The required length depends on the model:
#'   \describe{
#'     \item{\code{"constant"}}{Length 1: \code{c(k1)}}
#'     \item{\code{"constant_to"}}{Length 2: \code{c(k1, age_end)}}
#'     \item{\code{"constant_from"}}{Length 2: \code{c(k1, age_start)}}
#'     \item{\code{"constant_interval"}}{Length 3: \code{c(k1, age_start, age_end)}}
#'   }
#'   where \code{k1} is the transition rate (hazard coefficient, not probability).
#'
#' @param dx Numeric. The time step size in years. Determines the lesion
#'   acquisition probability (\code{k1 * dx}).
#'
#' @return A data.frame with the same structure as \code{state}, with updated
#'   values:
#'   \itemize{
#'     \item \code{lesion}: Updated to TRUE for agents who acquired a lesion.
#'     \item All other columns: Unchanged (age is not modified).
#'   }
#'
#' @examples
#' # Create a small test population
#' state <- data.frame(
#'   agent_id = 1:10,
#'   age = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
#'   lesion = rep(FALSE, 10)
#' )
#'
#' # Apply lesion acquisition with window [0, 6)
#' result <- apply_lesions(
#'   state = state,
#'   lesion_model = "constant_to",
#'   lesion_param = c(0.05, 6),  # k1 = 0.05, window ends at age 6
#'   dx = 1
#' )
#'
#' @seealso \code{\link{apply_mortality}} for the mortality module.
#'
#' @export
apply_lesions <- function(state,
                          lesion_model,
                          lesion_param,
                          dx) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (is.null(lesion_model)) {
    stop("lesion_model must be specified")
  }

  # ---------------------------------------------------------------------------
  # Extract parameters (model-specific)
  # ---------------------------------------------------------------------------

  if (lesion_model == "constant") {
    if (length(lesion_param) != 1) {
      stop("For constant model, lesion_param must be length 1: c(k1)")
    }
    k1 <- lesion_param[1]
    age_start <- -Inf
    age_end <- Inf

  } else if (lesion_model == "constant_to") {
    if (length(lesion_param) != 2) {
      stop("For constant_to model, lesion_param must be length 2: c(k1, age_end)")
    }
    k1 <- lesion_param[1]
    age_start <- -Inf
    age_end <- lesion_param[2]

  } else if (lesion_model == "constant_from") {
    if (length(lesion_param) != 2) {
      stop("For constant_from model, lesion_param must be length 2: c(k1, age_start)")
    }
    k1 <- lesion_param[1]
    age_start <- lesion_param[2]
    age_end <- Inf

  } else if (lesion_model == "constant_interval") {
    if (length(lesion_param) != 3) {
      stop("For constant_interval model, lesion_param must be length 3: c(k1, age_start, age_end)")
    }
    k1 <- lesion_param[1]
    age_start <- lesion_param[2]
    age_end <- lesion_param[3]

  } else {
    stop("lesion_model must be one of: 'constant', 'constant_to', 'constant_from', 'constant_interval'")
  }

  # ---------------------------------------------------------------------------
  # Copy state to avoid modifying input
  # ---------------------------------------------------------------------------

  result <- state

  # ---------------------------------------------------------------------------
  # Identify eligible agents
  # ---------------------------------------------------------------------------

  # Eligible if: alive, no lesion yet, and within age window
  eligible_idx <- which(
    !state$dead &
    !state$lesion &
    state$age >= age_start &
    state$age < age_end
  )

  # If no one is eligible, return early
  if (length(eligible_idx) == 0) {
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Calculate transition probability
  # ---------------------------------------------------------------------------

  p_lesion <- k1 * dx

  # ---------------------------------------------------------------------------
  # Determine who acquires lesions
  # ---------------------------------------------------------------------------

  u <- runif(length(eligible_idx))
  acquires <- u < p_lesion
  new_lesion_idx <- eligible_idx[acquires]

  # ---------------------------------------------------------------------------
  # Update state
  # ---------------------------------------------------------------------------

  result$lesion[new_lesion_idx] <- TRUE

  # Note: Age is NOT modified - that is the mortality module's responsibility

  # ---------------------------------------------------------------------------
  # Return updated state
  # ---------------------------------------------------------------------------

  return(result)
}

