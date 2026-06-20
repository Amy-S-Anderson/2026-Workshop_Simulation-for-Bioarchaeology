#' @title Apply Estimation Error to Agent Population
#'
#' @description Applies age estimation error to agents in the archaeological
#'   sample. Simulates the measurement error inherent in osteological age
#'   estimation methods, where uncertainty increases with age and systematic
#'   bias causes underestimation of older individuals.
#'
#' @details This function implements the final stage of archaeological
#'   filtration, adding realistic age estimation error to recovered skeletal
#'   remains. The error model has two components:
#'
#'   \enumerate{
#'     \item \strong{Random error}: Standard deviation increases with age.
#'       SD = 0 at age 0, rises linearly to \code{sd_at_20} at age 20, then
#'       continues increasing at \code{sd_per_year} for ages beyond 20.
#'     \item \strong{Systematic bias}: Older individuals are systematically
#'       underestimated. Bias = 0 up to \code{bias_start}, then increases
#'       linearly, reaching \code{bias_at_90} at age 90.
#'   }
#'
#'   Only agents who are in the sample (\code{in_sample = TRUE}) receive
#'   estimated ages. Agents not in the sample have \code{estimated_age = NA}.
#'
#'   The true age (\code{age}) is preserved; the estimated age is stored in
#'   a new column \code{estimated_age}.
#'
#' @param state A data.frame representing the current population state. Must
#'   contain the following columns:
#'   \describe{
#'     \item{agent_id}{Integer. Unique identifier for each agent.}
#'     \item{age}{Numeric. True age at death for dead agents.}
#'     \item{lesion}{Logical. TRUE if the agent has a lesion.}
#'     \item{dead}{Logical. TRUE if the agent is dead, FALSE if alive.}
#'     \item{was_deposited}{Logical. TRUE if the agent was deposited.}
#'     \item{in_sample}{Logical. TRUE if the agent is in the sample.}
#'   }
#'
#' @param error_model Character. The estimation error model to use. Currently
#'   only \code{"bespoke_increasing_sd"} is supported.
#'
#' @param error_param Numeric vector containing the error parameters. For the
#'   \code{"bespoke_increasing_sd"} model, this must be a length-4 vector:
#'   \describe{
#'     \item{sd_per_year}{(Element 1) Rate at which SD increases per year after age 20.
#'       Default calibration: 0.1375 gives SD ≈ 7.5 at age 60.}
#'     \item{sd_at_20}{(Element 2) SD of random error at age 20. Default: 2.}
#'     \item{bias_start}{(Element 3) Age at which systematic bias begins. Default: 60.}
#'     \item{bias_at_90}{(Element 4) Mean systematic error at age 90. Default: 25.}
#'   }
#'
#' @return A data.frame with the same structure as \code{state}, plus a new
#'   column:
#'   \itemize{
#'     \item \code{estimated_age}: Numeric. The estimated age with measurement
#'       error applied. NA for agents not in sample.
#'   }
#'   All original columns are preserved unchanged.
#'
#' @examples
#' # Create a test population of sampled individuals
#' state <- data.frame(
#'   agent_id = 1:10,
#'   age = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90),
#'   lesion = rep(FALSE, 10),
#'   dead = rep(TRUE, 10),
#'   was_deposited = rep(TRUE, 10),
#'   in_sample = rep(TRUE, 10)
#' )
#'
#' # Default error parameters
#' error_param <- c(0.1375, 2, 60, 25)
#'
#' # Apply estimation error
#' set.seed(42)
#' result <- apply_estimation_error(
#'   state = state,
#'   error_model = "bespoke_increasing_sd",
#'   error_param = error_param
#' )
#'
#' @seealso \code{\link{apply_preservation}} for the preservation filtering
#'   module that precedes estimation error.
#'
#' @export
apply_estimation_error <- function(state,
                                   error_model = "bespoke_increasing_sd",
                                   error_param = c(0.1375, 2, 60, 25)) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (is.null(error_model)) {
    stop("error_model must be specified")
  }

  # ---------------------------------------------------------------------------

  # Extract parameters (model-specific)
  # ---------------------------------------------------------------------------

  if (error_model == "bespoke_increasing_sd") {
    if (length(error_param) != 4) {
      stop("For bespoke_increasing_sd model, error_param must be length 4: ",
           "c(sd_per_year, sd_at_20, bias_start, bias_at_90)")
    }
    sd_per_year <- error_param[1]
    sd_at_20 <- error_param[2]
    bias_start <- error_param[3]
    bias_at_90 <- error_param[4]

  } else {
    stop("error_model must be 'bespoke_increasing_sd'. ",
         "Other models are not yet supported.")
  }

  # ---------------------------------------------------------------------------
  # Copy state to avoid modifying input
  # ---------------------------------------------------------------------------

  result <- state

  # ---------------------------------------------------------------------------
  # Initialize estimated_age column with NA
  # ---------------------------------------------------------------------------

  result$estimated_age <- rep(NA_real_, nrow(result))

  # ---------------------------------------------------------------------------
  # Identify eligible agents (in_sample = TRUE)
  # ---------------------------------------------------------------------------

  eligible_idx <- which(state$in_sample)

  # If no one is eligible, return early
  if (length(eligible_idx) == 0) {
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Apply estimation error (model-specific)
  # ---------------------------------------------------------------------------

  if (error_model == "bespoke_increasing_sd") {
    # Get true ages of eligible agents
    ages <- state$age[eligible_idx]

    # --- Random error: SD is piecewise linear in age ---
    # SD = 0 at age 0, rises linearly to sd_at_20 at age 20,
    # then continues rising at sd_per_year beyond age 20.
    sd_vec <- dplyr::case_when(
      ages <= 0  ~ 0,
      ages <= 20 ~ (ages / 20) * sd_at_20,
      TRUE       ~ sd_at_20 + (ages - 20) * sd_per_year
    )

    random_error <- stats::rnorm(length(ages), mean = 0, sd = sd_vec)

    # --- Systematic error: linear ramp starting at bias_start ---
    # Bias = 0 up to bias_start, then rises linearly.
    # Calibrated so bias = bias_at_90 at age 90.
    bias_slope <- bias_at_90 / (90 - bias_start)

    systematic_error <- dplyr::case_when(
      ages <= bias_start ~ 0,
      TRUE               ~ (ages - bias_start) * bias_slope
    )

    # Calculate estimated age
    estimated_ages <- ages + random_error - systematic_error

    # Store in result
    result$estimated_age[eligible_idx] <- estimated_ages
  }

  # ---------------------------------------------------------------------------
  # Return updated state
  # ---------------------------------------------------------------------------

  return(result)
}

