# --------------------------------------------------------------------------------------------------------------
# PERSEPHONE: Agent-Based Model for Bioarchaeology
# --------------------------------------------------------------------------------------------------------------
#
# This file contains the core simulation engine and its helper functions.
# Individual agents:
#   - Are born as a starting cohort of a specified size
#   - Face a stable annual risk of forming a skeletal lesion within a specified age range
#   - May experience different age-specific mortality depending on lesion status
#   - Age-specific mortality risks follow a Siler function


# --- Internal helpers (not exported) ---

#' Create the initial cohort data frame
#' @param cohort_size Number of agents in the starting cohort
#' @return A data frame with columns: agent_id, age, lesion, dead, in_sample
#' @keywords internal
create_cohort <- function(cohort_size) {
  data.frame(agent_id = 1:cohort_size,
             age = 0,
             lesion = 0,
             dead = FALSE,
             was_deposited = FALSE,
             in_sample = TRUE
  )
}

#' Age up all living agents to the current timestep
#' @param cohort Population data frame
#' @param k Current timestep (age to assign)
#' @return Updated cohort data frame
#' @keywords internal
age_cohort <- function(cohort, k) {
  Alive <- which(!cohort$dead)
  cohort$age[Alive] <- k
  cohort
}

#' Compute Siler hazard for a given age
#' @param k Age (integer timestep)
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric hazard value
#' @keywords internal
compute_siler_risk <- function(k, mortality_regime) {
  mortality_regime$a1 * exp(-mortality_regime$b1 * k) +
    mortality_regime$a2 +
    mortality_regime$a3 * exp(mortality_regime$b3 * k)
}

#' Vectorized lesion formation across all living agents
#'
#' Replaces the per-agent form_lesion() loop with a single vectorized draw.
#' Each living agent gets its own uniform draw; conditions are evaluated
#' element-wise so individual-level heterogeneity is preserved.
#'
#' @param cohort Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param formation_window_opens Age at which lesions can start forming
#' @param formation_window_closes Age at which lesions stop forming
#' @param lesion_formation_rate Annual probability of lesion formation
#' @return Updated cohort data frame
#' @keywords internal
form_lesions <- function(cohort, Alive,
                             formation_window_opens,
                             formation_window_closes,
                             lesion_formation_rate) {
  # One draw per living agent
  stress <- runif(length(Alive), 0, 1)
  
  ages   <- cohort$age[Alive]
  in_window <- ages >= formation_window_opens & ages <= formation_window_closes
  
  # An agent acquires a lesion if: in the formation window, stress roll passes,
  # AND it does not already have a lesion (pmax preserves previously acquired lesions)
  new_lesion <- as.integer(in_window & stress <= lesion_formation_rate)
  cohort$lesion[Alive] <- pmax(cohort$lesion[Alive], new_lesion)
  
  cohort
}

#' Vectorized mortality across all living agents
#'
#'
#' @param cohort Population data frame
#' @param Alive Integer vector of row indices for living agents
#' @param age_based_risk Baseline Siler mortality risk at current age (scalar)
#' @param mortality_risk_type One of "proportional", "time_decreasing", "time_increasing"
#' @param relative_mortality_risk Multiplier for lesion-bearing individuals
#' @return Updated cohort data frame
#' @keywords internal
apply_mortality <- function(cohort, Alive,
                                age_based_risk,
                                mortality_risk_type,
                                relative_mortality_risk) {
  death_dice <- runif(length(Alive), 0, 1)
  
  ages   <- cohort$age[Alive]
  lesion <- cohort$lesion[Alive]
  
  # Compute the effective death threshold for each agent.
  # case_when() evaluates conditions row-wise (element-wise on vectors), so
  # each agent follows exactly one branch — logically equivalent to the
  # original nested ifelse() but readable and easy to extend.
  threshold <- dplyr::case_when(
    lesion == 0 ~
      age_based_risk,
    
    lesion == 1 & mortality_risk_type == "proportional" ~
      age_based_risk * relative_mortality_risk,
    
    lesion == 1 & mortality_risk_type == "time_decreasing" ~
      age_based_risk * relative_mortality_risk / ((ages / 10) + relative_mortality_risk),
    
    lesion == 1 & mortality_risk_type == "time_increasing" ~
      age_based_risk * ((ages / 10) + relative_mortality_risk) / relative_mortality_risk,
    
    .default = age_based_risk   # fallback: treat as no-lesion baseline
  )
  
  cohort$dead[Alive] <- cohort$dead[Alive] | (death_dice < threshold)
  
  cohort
}

#' Record survivor snapshot for the current timestep
#' @param cohort Population data frame
#' @param k Current timestep
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal
record_survivors <- function(cohort, k) {
  n_alive  <- sum(!cohort$dead & cohort$age == k)
  n_lesion <- sum(!cohort$dead & cohort$lesion == 1 & cohort$age == k)
  data.frame(Age = k,
             Alive = n_alive,
             Lesion = n_lesion,
             Lesion_perc = ifelse(n_alive == 0, NA, round(n_lesion / n_alive * 100, 1)))
}

#' Finalize the cohort into a cemetery
#'
#' Ages remaining survivors one final step and marks all as dead.
#' @param cohort Population data frame
#' @param k Last completed timestep
#' @return Data frame with all agents marked dead
#' @keywords internal
finalize_cemetery <- function(cohort, k) {
  k <- k + 1
  Alive <- which(!cohort$dead)
  cohort$age[Alive] <- k
  cohort$dead <- TRUE
  cohort
}


# --- Main simulation function (exported) ---

#' Simulate a cemetery using the Persephone ABM
#'
#' Runs the Persephone agent-based model: a birth cohort ages through time,
#' facing annual risks of skeletal lesion formation and Siler-model mortality.
#' Individuals with lesions may experience modified mortality risk. The
#' simulation ends when fewer than 10 individuals remain alive.
#'
#' @param cohort_size Integer. Number of individuals in the starting cohort.
#' @param lesion_formation_rate Numeric. Annual probability of developing a
#'   lesion (between 0 and 1).
#' @param formation_window_opens Numeric. Age at which lesions can start
#'   forming. Default 0.
#' @param formation_window_closes Numeric. Age at which new lesions stop forming.
#' @param mortality_risk_type Character. How lesions modify mortality:
#'   "proportional", "time_decreasing", or "time_increasing".
#' @param relative_mortality_risk Numeric. Mortality multiplier for individuals
#'   with lesions. 1 = no effect, 2 = double risk. Default 1.
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{individual_outcomes}{Data frame of all individuals with age at
#'       death and lesion status.}
#'     \item{survivors}{Data frame of survivor counts and lesion prevalence
#'       at each age.}
#'   }
#'
#' @examples
#' result <- Simulate_Cemetery(
#'   cohort_size = 500,
#'   lesion_formation_rate = 0.10,
#'   formation_window_closes = 5,
#'   mortality_regime = CoaleDemenyWestF5
#' )
#'
#' @export
Simulate_Cemetery <- function(cohort_size,
                              lesion_formation_rate,
                              formation_window_opens = 0,
                              formation_window_closes,
                              mortality_risk_type = "proportional",
                              relative_mortality_risk = 1,
                              mortality_regime,
                              deposition_param = 0,
                              taphonomy_regime,
                              loss_strength = 'no_decay',
                              age_noise = FALSE) {
  cohort <- create_cohort(cohort_size)
  
  k <- 0  # Initialize time counter
  
  # Set up table for survivor output
  survivors <- vector("list", 100)
  
  # As long as more than 10 people are alive
  while (sum(!cohort$dead) >= 10) {
    k <- k + 1  # Increment time
    cohort <- age_cohort(cohort, k)
    Alive <- which(!cohort$dead)
    
    # Baseline Siler mortality risk for this age (scalar — same for all agents)
    age_based_risk <- compute_siler_risk(k, mortality_regime)
    
    # Vectorized lesion formation and mortality across all living agents.
    # Each agent receives its own independent draw; case_when() dispatches
    # each agent through the correct mortality branch element-wise.
    cohort <- form_lesions_vec(cohort, Alive,
                               formation_window_opens,
                               formation_window_closes,
                               lesion_formation_rate)
    
    cohort <- apply_mortality_vec(cohort, Alive,
                                  age_based_risk,
                                  mortality_risk_type,
                                  relative_mortality_risk)
    
    # Update summary log for survivors
    survivors[[k]] <- record_survivors(cohort, k)
  }
  
  # Once 10 or fewer people are left alive — they all enter the cemetery
  cohort <- finalize_cemetery(cohort, k)
  
  # Apply deposition bias (if any)
  cohort <- apply_deposition(cohort,
                             deposition_model = 'cutoff',
                             deposition_param = deposition_param,
                             dx = 1)
  cohort$in_sample <- cohort$was_deposited
  
  # Apply preservation bias (if any)
  if (loss_strength != 'no_decay') {
    a_siler <- c(taphonomy_regime$a1, taphonomy_regime$b1,
                 taphonomy_regime$a2, taphonomy_regime$a3,
                 taphonomy_regime$b3)
    b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
    cohort <- apply_preservation(cohort,
                                 preservation_model = 'siler',
                                 preservation_param = b_siler,
                                 dx = 1)
  }
  
  # Apply age misestimation (if any)
  if (age_noise) cohort <- apply_estimation_error(cohort)
  
  # Remove internal columns before returning
  cohort <- cohort %>% dplyr::select(-"dead")
  
  # Model output
  output <- list(individual_outcomes = cohort, survivors = do.call(rbind, survivors))
  
  return(output)
}