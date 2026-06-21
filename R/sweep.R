

#' Coale-Demeny West Female Level 5 mortality regime
#'
#' A mortality table used as the default regime in Persephone simulations.
#'
#' @name CoaleDemenyWestF5
#' @docType data
#' @keywords datasets
NULL
load('data/CoaleDemenyWestF5.rda')

#' @title Get default simulation parameters
#' @description Save this function, without any specified arguments, to an object, and that object will be a list of all the default parameter values for Simulate_Cemetery. 
#' @return A named list of default parameter values
#' @export
get_default_params <- function() {
  list(
    dx             = 1,
    max_years      = 100,
    pop0_size      = 1000,
    age_structured = FALSE,
    tfr            = NULL,
    mortality_regime = CoaleDemenyWestF5,
    
    lesion_formation_rate    = NULL,
    annual_exposure          = 0.1,
    lesion_formation_window  = c(0, 6),
    mortality_risk_type      = "proportional",
    lesion_related_hazard    = 1,
    
    gammafrailty_variance = NULL,

    exposure_causes_hazard     = FALSE,
    hazard_is_transient        = FALSE,
    lesion_requires_survival   = FALSE,
    exposure_hazard_multiplier = 1,
    
    deposition_param = 0,
    taphonomy_regime = NULL,
    loss_strength    = 'no_decay',
    age_noise        = FALSE,
    
    yearly_updates = FALSE
  )
}


#' Generate a directory name from a named list of scenario parameters
#'
#' Concatenates each parameter name and value as "param=value", separated by
#' underscores. Suitable for use as a folder name that is self-documenting.
#'
#' @param scenario_params Named list of parameters that vary across scenarios
#' @return Character string, e.g. "exposure_causes_hazard=TRUE_hazard_is_transient=FALSE_lesion_requires_survival=FALSE"
#' @export
scenario_dir_name <- function(scenario_params) {
  paste(
    mapply(function(val, nm) paste0(nm, "=", val),
           scenario_params, names(scenario_params)),
    collapse = "_"
  )
}


#' Run Simulate_Cemetery once and save all outputs to disk
#'
#' Creates a uniquely named subdirectory, saves the random seed, the full
#' parameter list as JSON, individual_outcomes as sim_cemetery.csv, and
#' annual_census as sim_census.csv.
#'
#' @param params Named list of parameters passed directly to Simulate_Cemetery
#' @param output_directory Full path to the directory for this single run
#' @param seed Optional integer seed. If NA a microsecond-based seed is chosen.
#' @return The simulation output list (individual_outcomes, annual_census),
#'   invisibly
#' @export
run_model_single <- function(params, output_directory, seed = NA) {
  
  if (dir.exists(output_directory)) {
    message("Directory already exists, skipping: ", output_directory)
    return(invisible(NULL))
  }
  dir.create(output_directory, recursive = TRUE)
  
  if (is.na(seed)) {
    seed <- as.numeric(format(Sys.time(), "%OS6")) * 1000000
  }
  set.seed(seed)
  write(seed, file = file.path(output_directory, "seed.csv"))
  
  # jsonlite can't serialise a data frame nested inside a list cleanly,
  # so store the mortality regime name as a string for the JSON record.
  params_for_json <- params
  params_for_json[["mortality_regime"]] <- params[["mortality_regime"]]$name
  write(jsonlite::toJSON(params_for_json, pretty = TRUE),
        file = file.path(output_directory, "params_used.json"))
  
  output <- do.call(Simulate_Cemetery, params)
  
  write.csv(output$individual_outcomes,
            file = file.path(output_directory, "sim_cemetery.csv"),
            row.names = FALSE)
  write.csv(output$annual_census,
            file = file.path(output_directory, "sim_census.csv"),
            row.names = FALSE)
  
  invisible(output)
}


#' Run repeated replicates of a single parameter set
#'
#' Each replicate is saved in its own rep=N subdirectory under
#' root_output_directory.
#'
#' @param params Named list of parameters for Simulate_Cemetery
#' @param root_output_directory Parent directory for this parameter set
#' @param numreps Integer number of replicates
#' @export
run_model_reps <- function(params, root_output_directory, numreps = 100) {
  for (rep in seq_len(numreps)) {
    output_directory <- file.path(root_output_directory, paste0("rep=", rep))
    run_model_single(params, output_directory)
  }
}


#' Run and save all scenario replicates, then return combined outputs
#'
#' For each scenario, merges the scenario parameters into base_params,
#' constructs a self-documenting directory name from the scenario parameters,
#' runs numreps replicates, then reads back and combines both output tables.
#'
#' Directory structure produced:
#'   root_output_directory/
#'     scenario=1_exposure_causes_hazard=FALSE_.../
#'       rep=1/
#'         sim_cemetery.csv
#'         sim_census.csv
#'         seed.csv
#'         params_used.json
#'       rep=2/
#'       ...
#'     scenario=2_.../
#'     ...
#'
#' @param base_params Named list of default parameters from get_default_params()
#' @param scenarios Named list of lists; each inner list holds the parameters
#'   that differ from base_params for that scenario
#' @param root_output_directory Root path under which all output is saved
#' @param numreps Integer number of replicates per scenario
#' @return A list with two combined data frames:
#'   \describe{
#'     \item{individual_outcomes}{All sim_cemetery rows, with columns
#'       scenario and rep appended}
#'     \item{annual_census}{All sim_census rows, with columns
#'       scenario and rep appended}
#'   }
#' @export
run_scenario_sweep <- function(base_params, scenarios,
                               root_output_directory,
                               numreps = 100) {
  
  all_cemetery <- vector("list", length(scenarios))
  all_census   <- vector("list", length(scenarios))
  
  for (i in seq_along(scenarios)) {
    scenario_name   <- names(scenarios)[i]
    scenario_params <- scenarios[[i]]
    
    # Merge scenario overrides into base params
    params <- base_params
    for (nm in names(scenario_params)) {
      params[[nm]] <- scenario_params[[nm]]
    }
    
    # Build a self-documenting directory name
    scenario_subdir <- paste0(
      "scenario=", scenario_name, "_",
      scenario_dir_name(scenario_params)
    )
    scenario_dir <- file.path(root_output_directory, scenario_subdir)
    
    message("Running scenario ", scenario_name, " (", numreps, " reps) ...")
    run_model_reps(params, scenario_dir, numreps = numreps)
    
    # Read back and combine replicates for this scenario
    rep_dirs <- list.dirs(scenario_dir, recursive = FALSE, full.names = TRUE)
    
    cemetery_reps <- lapply(seq_along(rep_dirs), function(j) {
      f <- file.path(rep_dirs[j], "sim_cemetery.csv")
      if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
      df <- readr::read_csv(f, show_col_types = FALSE)
      df$scenario <- scenario_name
      df$rep      <- j
      # Attach scenario parameter values as columns for easy filtering later
      for (nm in names(scenario_params)) df[[nm]] <- scenario_params[[nm]]
      df
    })
    
    census_reps <- lapply(seq_along(rep_dirs), function(j) {
      f <- file.path(rep_dirs[j], "sim_census.csv")
      if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
      df <- readr::read_csv(f, show_col_types = FALSE)
      df$scenario <- scenario_name
      df$rep      <- j
      for (nm in names(scenario_params)) df[[nm]] <- scenario_params[[nm]]
      df
    })
    
    all_cemetery[[i]] <- dplyr::bind_rows(Filter(Negate(is.null), cemetery_reps))
all_census[[i]]   <- dplyr::bind_rows(Filter(Negate(is.null), census_reps))
  }
  
  list(
    individual_outcomes = dplyr::bind_rows(all_cemetery),
    annual_census       = dplyr::bind_rows(all_census)
  )
}




### Function to read in the simulated data. This creates a nested list. 
read_scenario_sweep = function(root_output_directory, target_param, target_param_values, data_file) { # data_file is either "sim_cemetery.csv" or "sim_survivors.csv"
  sweep_results = list()
  
  # for each value in the parameter being swept
  for (param_value in target_param_values) {
    param_dir = file.path(root_output_directory, paste0(target_param,"=", param_value)) # identify the directory where the output for that scenario is stored
    rep_dirs = list.dirs(param_dir, recursive = FALSE, full.names = TRUE) # identify the file paths for all the iterations/runs of that scenario
    
    # extract the output .csv file from each run
    rep_outputs = lapply(rep_dirs, function(rep_dir) {
      output_file = file.path(rep_dir, data_file) 
      if (file.exists(output_file)) {
        return(read_csv(output_file))
      } else {
        print("file does not exist")
        return(NULL)
      }
    })
    
    # and save it as an entry in the sweep_results list, named with the target parameter value for this set of runs. 
    sweep_results[[paste0(target_param,"=",param_value)]] = rep_outputs
  }
  
  # return a list of lists: one for each value of the target parameter, containing a list of model output data frames from the output.csv files
  return(sweep_results)
}