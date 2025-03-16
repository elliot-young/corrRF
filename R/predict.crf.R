#' Predictions from a crf given newdata
#' @description Predictions from a fitted \code{crf} clustered random forest on newdata \code{newdata}.
#' @param object a fitted \code{crf} clustered random forest object fitted by \code{crf}.
#' @param newdata dataset on which predictions are to be performed.
#' @param sderr whether 'bootstrap of little bags' standard errors should be additionally outputted. Default is \code{FALSE}.
#' @param ... additional arguments
#' @return Fitted values, potentially alongside standard errors (see \code{sderr}).
#' @export
predict.crf <- function(object, newdata, sderr=FALSE, ...) {
  L <- length(object$forest)
  B <- length(object$forest[[1]])
  predictions <- lapply(object$forest, function(sublist) {
    lapply(sublist, function(x) {
      predict(x, newdata)
    })
  })
  prediction_mean_B <- lapply(predictions, function(Bs) Reduce("+", Bs)/B)
  mus <- colSums(do.call(rbind, prediction_mean_B))/L
  if (!sderr) {
    return(mus)
  } else {
    if (L<=10) warning("Warining: For bag of little bag bootstrap estimator of variance L should be large. If no variance estimator is required, set L = NULL (or 1)")
    vs <- rowSums(matrix((sapply(prediction_mean_B, identity) - mus)^2, nrow(newdata)) )/L
    return(list(fitted=mus, sderr=sqrt(vs)))
  }
}
