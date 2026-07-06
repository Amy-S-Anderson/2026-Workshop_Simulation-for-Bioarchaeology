#' Build a trauma-proportion lookup table and merge it onto a mu0 lookup
#'
#' Given a small set of user-specified control points describing what
#' proportion of all-cause mortality at a given age is attributable to
#' trauma (as opposed to illness), this function expands those control
#' points onto the same integer-age grid used by \code{defrail_siler()}'s
#' \code{mu0} lookup table, and returns the merged table for use in
#' \code{apply_mortality()}.
#'
#' Control points are linearly interpolated between specified ages. Ages
#' younger than the youngest control point, or older than the oldest, are
#' held flat at the nearest specified value (i.e., no extrapolation beyond
#' the given range). This means the analyst's assumptions about the oldest
#' and youngest ages are explicit rather than a modeling default silently
#' applied at the edges: e.g. if you supply a control point of
#' \code{p_trauma = 0.05} at age 70 and your model runs individuals to age
#' 100, ages 70-100 all get \code{p_trauma = 0.05} unless you add another
#' control point further out.
#'
#' @param mu0_lookup Data frame produced by \code{defrail_siler()}, with (at
#'   minimum) an integer \code{age} column.
#' @param control_points Data frame with columns \code{age} and
#'   \code{p_trauma}, giving the proportion of all-cause mortality due to
#'   trauma at each specified age. \code{p_trauma} values must fall 
#'   between 0 and 1 (inclusive). Ages need not be integers or cover the full range in
#'   \code{mu0_lookup}; they are interpolated/extended as described above.
#'   If \code{NULL} (the default), \code{p_trauma} is set to 0 at every age,
#'   i.e. all deaths are attributed to illness and no cause-of-death
#'   splitting occurs. 
#'
#' @return \code{mu0_lookup} with an added/overwritten \code{p_trauma}
#'   column, aligned to the same \code{age} grid.
#'
#' @examples
#' # Illness-dominant in infancy/childhood and old age; trauma-dominant
#' # from age 8 through middle adulthood:
#' control_points <- data.frame(
#'   age      = c(0,   4,   8,   15,  30,  50,  70,  90),
#'   p_trauma = c(0.05,0.10,0.40,0.55,0.55,0.35,0.15,0.05)
#' )
#'
#' @export
build_p_trauma_lookup <- function(mu0_lookup, control_points = NULL) {
  
  if (!"age" %in% names(mu0_lookup)) {
    stop("mu0_lookup must contain an 'age' column.")
  }
  
  ages <- mu0_lookup$age
  
  if (is.null(control_points)) {
    mu0_lookup$p_trauma <- 0
    return(mu0_lookup)
  }
  
  if (!all(c("age", "p_trauma") %in% names(control_points))) {
    stop("control_points must be a data frame with columns 'age' and 'p_trauma'.")
  }
  
  if (any(control_points$p_trauma < 0 | control_points$p_trauma > 1, na.rm = TRUE)) {
    stop("control_points$p_trauma values must all fall between 0 and 1.")
  }
  
  if (anyNA(control_points$age) || anyNA(control_points$p_trauma)) {
    stop("control_points must not contain NA values in 'age' or 'p_trauma'.")
  }
  
  # Sort by age and de-duplicate (keep the last value supplied for any
  # repeated age, so a user can overwrite a control point by re-specifying it)
  control_points <- control_points[order(control_points$age), ]
  control_points <- control_points[!duplicated(control_points$age, fromLast = TRUE), ]
  
  if (nrow(control_points) == 1) {
    # Single control point: flat p_trauma at every age
    p_trauma <- rep(control_points$p_trauma, length(ages))
  } else {
    # Linear interpolation between control points; flat extrapolation
    # (rule = 2) beyond the min/max control-point ages
    p_trauma <- stats::approx(
      x    = control_points$age,
      y    = control_points$p_trauma,
      xout = ages,
      rule = 2
    )$y
  }
  
  mu0_lookup$p_trauma <- p_trauma
  mu0_lookup
}
