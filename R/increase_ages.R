#' @title Increase all agent ages by one time step.
#'
#' @description Updates the 'age' column of the population data frame by adding 1 to all agent ages. Happy birthday!
#'
#' 
#' @param pop Population data frame
#' @return Updated pop data frame
#' @keywords internal
#' @export
age_pop <- function(pop) {
  pop$age <- pop$age + 1
  pop
}

# This function is called inside the main loop of Simulate_Cemetery().