#' Master data frame for all object to cache
#'
#' This is an internal object which will be updated on every new release
#' library(ExperimentHub)
#' eh <- query(ExperimentHub(localHub=FALSE), c("sesameData", "v1.13.1"))
#' data.frame(name=eh$title, eh=names(eh))
#'
#' @name df_master
#' @docType data
#' @return master sheet of sesameData objects
NULL

