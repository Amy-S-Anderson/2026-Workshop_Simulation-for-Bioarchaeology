# check_dx_resolution.R
#
# Standalone script to investigate how the time step (dx) affects the accuracy
# of persephone's forward simulation. Compares against demohaz's rejection
# sampler using two-sample Kolmogorov-Smirnov tests.
#
# Also tests both orderings of operations (lesions first vs mortality first)
# to verify correctness.
#
# Usage: source("check_dx_resolution.R") or run interactively

library(demohaz)
devtools::load_all()  # Load persephone from source

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

N <- 80000
dx_values <- c(0.01, 0.1, 0.5, 1.0)  # dx values to test
orders <- c("lesions_first", "mortality_first")  # operation orderings to test
dx_rej <- 0.001   # Rejection sampler grid spacing (reference)
xmax <- 120       # Maximum age
seed <- 12345     # Random seed for reproducibility

# -----------------------------------------------------------------------------
# Model parameters (matching test-usher3-functional.R)
# -----------------------------------------------------------------------------

b_siler_test <- c(.175, 1.40, .368 * .01,
                  log(.917 * .1 / (.075 * .001)) / (.917 * .1),
                  .917 * .1)
th0 <- c(2e-2, 1.2, b_siler_test)

k1 <- th0[1]
k2 <- th0[2]
b_siler <- th0[3:7]
mortality_param <- c(k2, b_siler)

# -----------------------------------------------------------------------------
# Forward simulation helper (from test-usher3-functional.R)
# -----------------------------------------------------------------------------

simulate_persephone_forward <- function(N, k1, mortality_param, dx, xmax,
                                         lesion_model = "constant",
                                         x_cut = Inf,
                                         order = "lesions_first") {
  if (lesion_model == "constant") {
    lesion_param <- c(k1)
  } else if (lesion_model == "constant_to") {
    lesion_param <- c(k1, x_cut)
  } else {
    stop("Unsupported lesion_model: ", lesion_model)
  }

  if (!order %in% c("lesions_first", "mortality_first")) {
    stop("order must be 'lesions_first' or 'mortality_first'")
  }

  state <- data.frame(
    agent_id = 1:N,
    age = rep(0, N),
    lesion = rep(FALSE, N),
    dead = rep(FALSE, N),
    in_sample = rep(TRUE, N)
  )

  nsteps <- ceiling(xmax / dx)

  for (step in seq_len(nsteps)) {
    if (all(state$dead)) break

    alive_before <- !state$dead

    if (order == "lesions_first") {
      # Lesions first, then mortality
      state <- apply_lesions(
        state = state,
        lesion_model = lesion_model,
        lesion_param = lesion_param,
        dx = dx
      )

      state <- apply_mortality_usher3(
        state = state,
        mortality_model = "usher3",
        mortality_param = mortality_param,
        dx = dx
      )
    } else {
      # Mortality first, then lesions
      state <- apply_mortality_usher3(
        state = state,
        mortality_model = "usher3",
        mortality_param = mortality_param,
        dx = dx
      )

      state <- apply_lesions(
        state = state,
        lesion_model = lesion_model,
        lesion_param = lesion_param,
        dx = dx
      )
    }

    died_this_step <- alive_before & state$dead
    n_died <- sum(died_this_step)
    if (n_died > 0) {
      state$age[died_this_step] <- state$age[died_this_step] + runif(n_died, 0, dx)
    }
  }

  list(
    x = state$age,
    ill = state$lesion
  )
}

# -----------------------------------------------------------------------------
# Run comparison for each dx value and order
# -----------------------------------------------------------------------------

cat("=============================================================\n")
cat("dx Resolution Check for Persephone Forward Simulation\n")
cat("=============================================================\n")
cat(sprintf("N = %d samples per method\n", N))
cat(sprintf("dx values to test: %s\n", paste(dx_values, collapse = ", ")))
cat(sprintf("Orders to test: %s\n", paste(orders, collapse = ", ")))
cat(sprintf("Reference (demohaz) dx_rej = %g\n", dx_rej))
cat("=============================================================\n\n")

set.seed(seed)

# Generate reference sample once (demohaz rejection sampler)
cat("Generating reference sample (demohaz rejection sampler)...\n")
demo <- demohaz::sample_usher3(N, th0, dx_rej, xmax, x_cut = Inf)
demo_well <- demo$x[!demo$ill]
demo_ill <- demo$x[demo$ill]
prop_demo <- mean(demo$ill)
cat(sprintf("  Reference: %d well, %d ill (%.2f%% lesioned)\n\n",
            length(demo_well), length(demo_ill), 100 * prop_demo))

# Store results
results <- data.frame(
  dx = numeric(),
  order = character(),
  ks_p_well = numeric(),
  ks_p_ill = numeric(),
  prop_perse = numeric(),
  prop_diff = numeric(),
  prop_z = numeric(),
  pass_well = logical(),
  pass_ill = logical(),
  pass_prop = logical(),
  pass_all = logical(),
  stringsAsFactors = FALSE
)

for (dx in dx_values) {
  for (ord in orders) {
    cat(sprintf("Testing dx = %g, order = %s ...\n", dx, ord))

    perse <- simulate_persephone_forward(
      N = N,
      k1 = k1,
      mortality_param = mortality_param,
      dx = dx,
      xmax = xmax,
      lesion_model = "constant",
      x_cut = Inf,
      order = ord
    )

    perse_well <- perse$x[!perse$ill]
    perse_ill <- perse$x[perse$ill]

    ks_well <- ks.test(perse_well, demo_well)
    ks_ill <- ks.test(perse_ill, demo_ill)

    prop_perse <- mean(perse$ill)
    prop_diff <- abs(prop_perse - prop_demo)

    # SE for difference of two proportions (pooled estimate)
    p_pool <- (prop_perse + prop_demo) / 2
    se_diff <- sqrt(2 * p_pool * (1 - p_pool) / N)
    prop_z <- prop_diff / se_diff

    pass_well <- ks_well$p.value > 0.05
    pass_ill <- ks_ill$p.value > 0.05
    pass_prop <- prop_diff < 0.03
    pass_all <- pass_well && pass_ill && pass_prop

    results <- rbind(results, data.frame(
      dx = dx,
      order = ord,
      ks_p_well = ks_well$p.value,
      ks_p_ill = ks_ill$p.value,
      prop_perse = prop_perse,
      prop_diff = prop_diff,
      prop_z = prop_z,
      pass_well = pass_well,
      pass_ill = pass_ill,
      pass_prop = pass_prop,
      pass_all = pass_all,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  Persephone: %d well, %d ill (%.2f%% lesioned)\n",
                length(perse_well), length(perse_ill), 100 * prop_perse))
    cat(sprintf("  KS p-value (well): %.4f %s\n",
                ks_well$p.value, ifelse(pass_well, "PASS", "FAIL")))
    cat(sprintf("  KS p-value (ill):  %.4f %s\n",
                ks_ill$p.value, ifelse(pass_ill, "PASS", "FAIL")))
    cat(sprintf("  Prop diff:         %.4f (z=%.2f) %s\n",
                prop_diff, prop_z, ifelse(pass_prop, "PASS", "FAIL")))
    cat(sprintf("  Overall: %s\n\n", ifelse(pass_all, "PASS", "FAIL")))
  }
}

# -----------------------------------------------------------------------------
# Full Results Summary
# -----------------------------------------------------------------------------

cat("=============================================================\n")
cat("Full Results\n")
cat("=============================================================\n")
print(results[, c("dx", "order", "ks_p_well", "ks_p_ill", "prop_diff", "prop_z", "pass_all")])
cat("\n")

# Find largest dx that passes for each order
cat("Largest dx that passes all checks:\n")
for (ord in orders) {
  ord_results <- results[results$order == ord, ]
  passing_dx <- ord_results$dx[ord_results$pass_all]
  if (length(passing_dx) > 0) {
    cat(sprintf("  %s: dx = %g\n", ord, max(passing_dx)))
  } else {
    cat(sprintf("  %s: none\n", ord))
  }
}

# -----------------------------------------------------------------------------
# Lesion Proportion Summary Table
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("Lesion Proportion Summary\n")
cat("=============================================================\n")
cat(sprintf("Reference (demohaz): %.4f (%.2f%%)\n", prop_demo, 100 * prop_demo))

se_expected <- sqrt(2 * prop_demo * (1 - prop_demo) / N)
cat(sprintf("Expected SE for difference: %.5f (%.3f%%)\n\n", se_expected, 100 * se_expected))

# Create clean proportion table
prop_table <- data.frame(
  dx = results$dx,
  order = results$order,
  prop = sprintf("%.4f", results$prop_perse),
  diff = sprintf("%+.4f", results$prop_perse - prop_demo),
  z = sprintf("%+.2f", sign(results$prop_perse - prop_demo) * results$prop_z),
  viable = ifelse(results$prop_z < 2, "yes", "NO")
)

print(prop_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# Order Effect Analysis
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("Order Effect Analysis\n")
cat("=============================================================\n")
cat("Comparing lesions_first vs mortality_first at each dx:\n\n")

for (dx in dx_values) {
  lf <- results[results$dx == dx & results$order == "lesions_first", ]
  mf <- results[results$dx == dx & results$order == "mortality_first", ]

  diff_prop <- lf$prop_perse - mf$prop_perse
  # SE for difference between two persephone runs (independent)
  se_order <- sqrt(2 * prop_demo * (1 - prop_demo) / N)
  z_order <- diff_prop / se_order

  cat(sprintf("dx = %g:\n", dx))
  cat(sprintf("  lesions_first:   %.4f (%.2f%%)\n", lf$prop_perse, 100 * lf$prop_perse))
  cat(sprintf("  mortality_first: %.4f (%.2f%%)\n", mf$prop_perse, 100 * mf$prop_perse))
  cat(sprintf("  difference:      %+.4f (z = %+.2f) %s\n\n",
              diff_prop, z_order,
              ifelse(abs(z_order) < 2, "", "<-- ORDER MATTERS")))
}
