# check_quantization_error.R
#
# Validates the one-step quantization error formula (eq. 9 in
# quantization_derivations.tex) by direct simulation.
#
# For a grid of (age, starting lesion proportion) pairs, the script:
#   1. Constructs a synthetic cohort with exactly that proportion lesioned.
#   2. Runs a single simulation step (lesions first, then mortality) with dx=1.
#   3. Repeats many times and averages the resulting proportion.
#   4. Compares the mean observed error (relative to the exact continuous-time
#      solution) against the closed-form predicted quantization error.
#
# This tests each grid point in isolation — no trajectory coupling between
# steps — so the comparison is a direct check of the one-step formula.
#
# Output: inst/scripts/quantization_error_dx1.png
#
# Usage: source("inst/scripts/check_quantization_error.R")

library(demohaz)
devtools::load_all()

# --- Configuration -----------------------------------------------------------

N      <- 10000  # agents per cohort per replication
n_reps <- 1000   # replications per grid point
h      <- 1      # step width (dx = 1)
seed   <- 12345

# Ages at which to evaluate (determines baseline Siler hazard mu).
# Range 0-80 spans the Siler U-shape; beyond ~80, mu*h approaches 1
# and the first-order approximation underlying the quantized update
# breaks down.
ages <- seq(0, 80, by = 10)

# Starting lesion proportions to test.  p = 1 is omitted because the
# error is trivially zero when everyone is already lesioned.
p_starts <- seq(0, 0.9, by = 0.05)

# Model parameters (matching check_dx_resolution.R / test-usher3-functional.R)
b_siler_test <- c(.175, 1.40, .368 * .01,
                  log(.917 * .1 / (.075 * .001)) / (.917 * .1),
                  .917 * .1)
th0 <- c(2e-2, 1.2, b_siler_test)

k1 <- th0[1]       # well-to-lesioned transition rate
k2 <- th0[2]       # mortality multiplier for lesioned individuals
b_siler <- th0[3:7] # Siler hazard parameters (demohaz parameterization)
mortality_param <- c(k2, b_siler)

# --- Theoretical helpers (from quantization_derivations.tex) -----------------

#' Auxiliary function psi(h, delta).
#'
#' Defined in the derivation as (1 - exp(-delta*h)) / delta when delta != 0,
#' or h when delta = 0.
#'
#' @param h     Numeric scalar. Step width.
#' @param delta Numeric vector. delta = k1 + (1 - k2) * mu.
#' @return Numeric vector, same length as delta.
psi <- function(h, delta) {
  ifelse(abs(delta) < 1e-12, h, (1 - exp(-delta * h)) / delta)
}

#' Exact continuous-time proportion lesioned among survivors after one step.
#'
#' Equation (5) in quantization_derivations.tex.
#'
#' @param p  Numeric. Starting proportion lesioned (in [0, 1]).
#' @param h  Numeric scalar. Step width.
#' @param k1 Numeric scalar. Well-to-lesioned transition rate.
#' @param k2 Numeric scalar. Mortality multiplier for lesioned individuals.
#' @param mu Numeric. Baseline (Siler) mortality hazard at the current age.
#' @return Numeric. Exact proportion lesioned among survivors at age + h.
p_exact <- function(p, h, k1, k2, mu) {
  delta <- k1 + (1 - k2) * mu
  ps <- psi(h, delta)
  num <- p + (1 - p) * k1 * ps
  den <- p + (1 - p) * (exp(-delta * h) + k1 * ps)
  num / den
}

#' Quantized (linearized) proportion lesioned among survivors after one step.
#'
#' Equation (8) in quantization_derivations.tex.  This is the expected outcome
#' of the split first-order update used in the forward simulation (lesions
#' first, then mortality).
#'
#' @param p  Numeric. Starting proportion lesioned (in [0, 1]).
#' @param h  Numeric scalar. Step width.
#' @param k1 Numeric scalar. Well-to-lesioned transition rate.
#' @param k2 Numeric scalar. Mortality multiplier for lesioned individuals.
#' @param mu Numeric. Baseline (Siler) mortality hazard at the current age.
#' @return Numeric. Quantized proportion lesioned among survivors at age + h.
p_quantized <- function(p, h, k1, k2, mu) {
  num <- (p + (1 - p) * k1 * h) * (1 - k2 * mu * h)
  den <- 1 - mu * ((1 - p) + p * k2) * h +
    (1 - p) * k1 * mu * (1 - k2) * h^2
  num / den
}

# --- Single-step simulation --------------------------------------------------

#' Run one step of the forward simulation on a synthetic cohort.
#'
#' Constructs a cohort of N agents at the given age, with exactly
#' floor(p * N) agents already lesioned, then applies one step of the
#' simulation (lesions first, mortality second) with step width h.
#'
#' @param N              Integer. Number of agents in the cohort.
#' @param age            Numeric scalar. Current age of all agents.
#' @param p              Numeric scalar. Starting proportion lesioned.
#' @param k1             Numeric scalar. Well-to-lesioned transition rate.
#' @param mortality_param Numeric vector of length 6: c(k2, b_siler[1:5]).
#' @param h              Numeric scalar. Step width.
#' @return Numeric scalar. Proportion lesioned among survivors after the step,
#'         or NA if no agents survive.
run_one_step <- function(N, age, p, k1, mortality_param, h) {
  n_lesioned <- floor(p * N)

  state <- data.frame(
    agent_id  = seq_len(N),
    age       = rep(age, N),
    lesion    = c(rep(TRUE, n_lesioned), rep(FALSE, N - n_lesioned)),
    dead      = rep(FALSE, N),
    in_sample = rep(TRUE, N)
  )

  # Lesions first, then mortality — same ordering as the forward simulation
  state <- apply_lesions(
    state        = state,
    lesion_model = "constant",
    lesion_param = k1,
    dx           = h
  )
  state <- apply_mortality_usher3(
    state           = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx              = h
  )

  alive <- !state$dead
  if (any(alive)) {
    sum(state$lesion[alive]) / sum(alive)
  } else {
    NA_real_
  }
}

# --- Sweep over the (age, p) grid --------------------------------------------

set.seed(seed)

# Build the grid: one row per (age, p_start) combination
grid <- expand.grid(age = ages, p_start = p_starts)
n_grid <- nrow(grid)

# Baseline Siler hazard at each grid age
grid$mu <- demohaz::hsiler(grid$age, b_siler)

# Theoretical predictions using the exact formulas
grid$pred_exact <- p_exact(grid$p_start, h, k1, k2, grid$mu)
grid$pred_quant <- p_quantized(grid$p_start, h, k1, k2, grid$mu)
grid$pred_error <- grid$pred_quant - grid$pred_exact

# Run replications and accumulate the observed proportion after one step
grid$mean_p_end <- NA_real_

for (i in seq_len(n_grid)) {
  cat(sprintf("\rGrid point %d / %d (age = %g, p = %.2f)",
              i, n_grid, grid$age[i], grid$p_start[i]))

  p_end_vec <- replicate(n_reps, run_one_step(
    N = N, age = grid$age[i], p = grid$p_start[i],
    k1 = k1, mortality_param = mortality_param, h = h
  ))

  # Drop NAs (all-dead cohorts — rare at these ages and cohort size)
  grid$mean_p_end[i] <- mean(p_end_vec, na.rm = TRUE)
}
cat("\n")

# Observed quantization error: mean simulated outcome minus exact prediction
grid$obs_error <- grid$mean_p_end - grid$pred_exact

# --- Scatter plot: predicted vs observed quantization error ------------------

#' Map ages to a blue-to-red colour ramp.
#'
#' @param ages    Numeric vector of ages.
#' @param max_age Numeric scalar. Age that maps to the red end of the ramp.
#' @return Character vector of hex colours, same length as ages.
age_to_col <- function(ages, max_age) {
  pal <- colorRampPalette(c("steelblue", "tomato"))
  idx <- pmax(1, round((ages / max_age) * 100))
  pal(100)[idx]
}

png("inst/scripts/quantization_error_dx1.png",
    width = 700, height = 700, res = 120)

lim <- range(c(grid$pred_error, grid$obs_error))
lim <- lim + diff(lim) * c(-0.05, 0.05)

plot(grid$pred_error, grid$obs_error,
     xlab = "Predicted quantization error (eq. 9)",
     ylab = paste0("Mean observed error (", n_reps, " reps \u00d7 N = ",
                    formatC(N, big.mark = ","), ")"),
     main = expression("Predicted vs observed quantization error," ~ dx == 1),
     xlim = lim, ylim = lim, asp = 1,
     pch = 16, col = age_to_col(grid$age, max(ages)), cex = 0.9)
abline(0, 1, lty = 2, col = "grey40")

n_ticks <- 5
tick_ages <- round(seq(0, max(ages), length.out = n_ticks))
legend("topleft",
       legend = paste("age", tick_ages),
       pch = 16,
       col = age_to_col(tick_ages, max(ages)),
       cex = 0.7,
       title = "Age")

dev.off()

cat("Plot saved to inst/scripts/quantization_error_dx1.png\n")
