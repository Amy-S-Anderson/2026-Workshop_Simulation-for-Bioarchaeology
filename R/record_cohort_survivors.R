
#' Helper function: Annual Census
#' 
#' Record survivor snapshot for the current timestep. 
#' @param pop Population data frame
#' @param current_time Current timestep
#' @param model_lesions Logical
#' @return One-row data frame with Age, Alive, Lesion, Lesion_perc
#' @keywords internal

record_cohort_survivors <- function(pop, current_time, model_lesions) {
  
  survivors <- data.frame(Time = current_time,
                          Age = unique(pop$age), 
                          Alive = nrow(pop))
  if(model_lesions){
    survivors$Lesion = sum(pop$lesion, na.rm = TRUE)
    survivors$Lesion_perc = round((sum(pop$lesion, na.rm = TRUE) / nrow(pop)) * 100, 1)
  }
  survivors
}

