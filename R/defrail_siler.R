#' Defrail a Siler mortality regime
#'
#' Given a Siler mortality regime whose parameters were fit to observed
#' (population-aggregate) mortality data, compute the underlying individual-
#' level baseline hazard mu_0(a) that, when mixed over a Gamma frailty
#' distribution with mean 1 and variance frailty_variance, reproduces the
#' observed hazard. Returns a lookup table of mu_0 by integer age.
#'
#' The Gamma-frailty mixture relationship is:
#'   mu_bar(a) = mu_0(a) / (1 + sigma^2 * H_0(a))
#' which implies the ODE:
#'   dH_0/da = mu_bar(a) * (1 + sigma^2 * H_0(a)),  H_0(0) = 0
#' solved numerically via RK4, then mu_0(a) = dH_0/da.
#'
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @param frailty_variance Numeric scalar, variance of the Gamma frailty
#'   distribution (mean is fixed at 1). Set to 0 to recover the observed
#'   hazard unchanged.
#' @param max_age Integer, oldest age to compute (default 120)
#' @param step Numeric, integration step size in years (default 0.01; finer
#'   steps improve accuracy at the cost of speed)
#'
#' @return A data frame with columns:
#'   \item{age}{Integer ages 0:max_age}
#'   \item{mu0}{Defrailted baseline hazard at each age}
#' @export
defrail_siler <- function(mortality_regime,
                          frailty_variance,
                          max_age = 120,
                          step    = 0.01) {
  
  if (frailty_variance < 0) {
    stop("frailty_variance must be non-negative")
  }
  
  # If frailty_variance == 0, the defrailed hazard is just the observed hazard
  if (frailty_variance == 0) {
    ages <- 0:max_age
    return(data.frame(
      age = ages,
      mu0 = compute_siler_risk(ages, mortality_regime)
    ))
  }
  
  s2 <- frailty_variance
  
  # RHS of the ODE: dH0/da = mu_bar(a) * (1 + s2 * H0)
  dH0_da <- function(a, H0) {
    mu_bar <- compute_siler_risk(a, mortality_regime)
    mu_bar * (1 + s2 * H0)
  }
  
  # Integrate via RK4 over a fine grid, then read off integer ages
  fine_ages <- seq(0, max_age, by = step)
  n         <- length(fine_ages)
  H0_fine   <- numeric(n)   # cumulative baseline hazard on fine grid
  mu0_fine  <- numeric(n)   # defrailed hazard on fine grid
  
  H0_fine[1]  <- 0
  mu0_fine[1] <- dH0_da(0, 0)
  
  for (i in seq_len(n - 1)) {
    a  <- fine_ages[i]
    H  <- H0_fine[i]
    
    k1 <- dH0_da(a,              H)
    k2 <- dH0_da(a + step / 2,   H + step * k1 / 2)
    k3 <- dH0_da(a + step / 2,   H + step * k2 / 2)
    k4 <- dH0_da(a + step,       H + step * k3)
    
    H0_fine[i + 1]  <- H + (step / 6) * (k1 + 2*k2 + 2*k3 + k4)
    mu0_fine[i + 1] <- dH0_da(a + step, H0_fine[i + 1])
  }
  
  # Sample mu0 at integer ages for the lookup table
  integer_indices <- round(0:max_age / step) + 1
  integer_indices <- pmin(integer_indices, n)   # guard against rounding past end
  
  data.frame(
    age = 0:max_age,
    mu0 = mu0_fine[integer_indices]
  )
}