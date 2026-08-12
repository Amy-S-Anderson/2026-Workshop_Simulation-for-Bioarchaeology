

#' List built-in taphonomy regimes (preservation loss)
#'
#' Scans the persephone package namespace and returns every object tagged
#' with class "persephone_taphonomy_regime", regardless
#' of what other classes it also carries (like "data.frame"). This means
#' any new regime added to the package with the right class tag will show
#' up here automatically — nothing needs to be manually registered.
#' Note that the taphonomy regimes are also Siler function parameters, since preservation loss, like mortality hazard, tends to be greatest at the youngest and oldest ages. 
#' 
#' 
#' @return A named list of Siler taphonomy regime objects (classed data
#'   frames), where the names match the object names as they exist in the
#'   package (e.g. "strong_decay").
#'
#' @export
list_taphonomy_regimes <- function() {
  
  # data(package = "persephone") asks R directly "what datasets does this
  # installed package provide?" — this works regardless of whether the
  # package has been attached via library(), because it reads the
  # package's data index rather than scanning a live environment.
  dataset_names <- data(package = "persephone")$results[, "Item"]
  
  # get() with envir/pos as an empty environment forces each dataset to
  # actually load (triggering its lazy-load promise) and returns it.
  all_objects <- mget(dataset_names, envir = as.environment("package:persephone"))
  
  is_regime <- vapply(
    all_objects,
    function(x) inherits(x, "persephone_taphonomy_regime"),
    FUN.VALUE = logical(1)
  )
  
  all_objects[is_regime]
}
