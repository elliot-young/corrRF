#' Summary for a crf fitted object
#' @description Summary of a fitted \code{crf} clustered random forest object fitted by \code{crf}.
#' @param object a fitted \code{crf} clustered random forest object fitted by \code{crf}.
#' @param ... additional arguments
#' @return Prints summary output for \code{crf} object
#' @export
summary.crf <- function(object, ...) {
  L <- length(object$forest)
  B <- length(object$forest[[1]])
  cat("crf: clustered random forest \n\n")
  cat("Call:\n", deparse(object$call), "\n\n")
  if (L==1) {
    cat("Number of trees:                 ", B, "\n")
  } else {
    cat("Number of bags for little bags:  ", L, "\n")
    cat("Number of little bags:           ", B, "\n")
  }
  cat("Number of observations:          ", object$num.samples, "\n")
  cat("Number of clusters/ groups:      ", object$num.groups, "\n")
  cat("Number of covariates:            ", object$num.of.covariates, "\n")
  cat("Subsampling rate:                ", object$beta, "\n")
  cat("Honesty type:                    ", ifelse(object$honesty, "honest", "dishonest"), "\n")
  cat("Working weight structure:        ", ifelse(object$correlation=="equicorr", "equicorrelated", "autoregressive AR(1)"), "\n")
  cat("Weight optimisation:             ", ifelse(object$honesty, "honest", "dishonest"), "\n")
}
