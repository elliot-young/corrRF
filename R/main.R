#' Clustered random forest fitting
#'
#' @param formula an object of class `formula` describing the model to fit.
#' @param data training dataset for fitting the CRF. Note that group ID must be given by the column \code{id}.
#' @param B the total number of trees (or trees per little bag if \eqn{L\neq}`NULL`). Default is 500.
#' @param L the total number of little bags if providing a bootstrap of little bags estimate for inference. To not include set \eqn{L=}`NULL`. Default is `NULL`.
#' @param beta the subsampling rate. Default is \eqn{beta=0.9}.
#' @param weight_optimiser the method used to construct weights. Options are `Pointwise variance`, `Training MSE` or `Test MSE`. Default is `Training MSE`.
#' @param correlation the weight structure implemented. Currently supported options are `ar1` and `equicorr`. Default is `equicorr`.
#' @param maxdepth the maximum depth of the decision tree fitting. Default is 30.
#' @param minbucket the minbucket of the decision tree fitting. Default is 10.
#' @param cp the complexity paramter for decision tree fitting. Default is 0.
#' @param x0 the covariate point to optimise weights towards if `weightoptimiser` set to `Pointwise variance`.
#' @param test_data the test dataset to optimise weights towards if `weightoptimiser` set to `Test MSE`.
#' @param fixrho fixes a pre-specified weight structure, given by the relevant `ar1` or `equicorr` parameter. Default is `FALSE` (optimise weights).
#' @param honesty whether honest or dishonest trees to be fit. Default is `TRUE`.
#'
#' @return A clustered random forest fitted object
#' @importFrom rpart rpart predict
#' @importFrom stats setNames terms
#' @export
crf <- function(formula, data, B=500, L=100, beta=0.9, weight_optimiser="Training MSE", correlation="equicorr", maxdepth=30, minbucket=10, cp=0, x0=NULL, test_data=NULL, fixrho=FALSE, honesty=TRUE) {

  if (!"id" %in% colnames(data)) stop("data requires an 'id' column that indicates cluster/ grouping structure")
  if (!is.null(L)) { if ((L<20) && (L>1)) warning("Warining: For bag of little bag bootstrap estimator of variance L should be large. If no variance estimator is required, set L = NULL (or 1)") }
  one_or_two = ifelse(is.null(L) || L == 1, 1, 2)
  L = ifelse(is.null(L) || L == 1, 1, L)
  optimise_randeff <- switch(paste(weight_optimiser, correlation),
                             "Pointwise variance equicorr" = optimise_randeff_ptwvar_equicorr,
                             "Pointwise variance ar1" = optimise_randeff_ptwvar_ar1,
                             "Training MSE equicorr" = optimise_randeff_trainmse_equicorr,
                             "Training MSE ar1" = optimise_randeff_trainmse_ar1,
                             "Test MSE equicorr" = optimise_randeff_testmse_equicorr,
                             "Test MSE ar1" = optimise_randeff_testmse_ar1,
                             stop("Error: Invalid weight optimiser and/or correlation! Weight optimiser must be either Pointwise variance, Training MSE or Test MSE. Correlation must be either equicorr or ar1.")
  )
  XWX_XWY_calc <- switch(correlation,
                         "equicorr" = XWX_XWY_equicorr_cpp,
                         "ar1" = XWX_XWY_ar1_cpp,
                         stop("Error: Invalid correlation! Correlation must be either equicorr or ar1.")
  )
  colnames(data)[colnames(data)==all.vars(formula)[1]] <- "y" # Set response variable to 'y'

  if (weight_optimiser == "Pointwise variance") {
    covariates <- all.vars(formula)[-1]
    if (length(x0) != length(covariates)) stop("Error: x0 of different length to number of covariates")
    x0dataset <- if (mode(x0) == "numeric") data.frame(t(x0), colnames = covariates) else x0
    if (!all(sort(colnames(x0dataset)) == sort(covariates))) {
      stop("Error: The datapoint x0 should be a dataframe with one row, with columns the covariates given in formula")
    }
  }

  forests <- list()
  starttime <- proc.time(); printedtime <- 1
  for (l in seq_len(L)) {
    set.seed(l)
    unique_ids <- unique(data$id)
    half_data <- data[data$id %in% sample(unique_ids, ceiling(length(unique_ids)/one_or_two)),]
    forest <- list()
    full_ids <- half_data$id
    unique_ids <- unique(full_ids)
    I <- length(unique_ids)
    s <- ceiling(I^beta)
    for (b in 1:B) {
      # Timings
      timeinmins <- (proc.time() - starttime)[[3]] / 60
      timeincriments <- ifelse(timeinmins > 1, max(0.1*(B*L/((l-1)*B+b))*timeinmins, 1), 0)
      if (timeinmins > printedtime) {
        print(paste0("Time Elapsed: ", floor(timeinmins), " minutes; Estimated Time Remaining: ", ceiling(timeinmins * (B * L / ((l - 1) * B + b) - 1)), " minutes."))
        printedtime <- printedtime + timeincriments
      }

      # Split data (honesty)
      set.seed((l-1)*B+b)
      ids_selected <- sample(unique_ids, s)
      rdf <- half_data[full_ids %in% ids_selected,]
      unique_ids_3partition <- split(ids_selected, cut(seq_along(ids_selected), breaks = 3, labels = FALSE))
      rdf_split <- rdf[rdf$id %in% unique_ids_3partition[[1]],]
      rdf_evalfix <- rdf[rdf$id %in% unique_ids_3partition[[2]],]
      rdf_evalrand <- rdf[rdf$id %in% unique_ids_3partition[[3]],]
      if (!honesty) rdf_split <- rdf_evalfix <- rdf_evalrand

      # Create initial (unweighted) tree using rpart
      initial_tree <- rpart(formula, rdf_split, cp=cp, maxdepth=maxdepth, minbucket=minbucket)
      # Leaf encoder
      leaves <- which(initial_tree$frame$var == "<leaf>")
      num_leaves <- length(leaves)
      initial_tree_where <- initial_tree$where
      initial_tree_where_1ton_indexed <- as.integer(initial_tree$where)
      lookup_preds <- setNames(as.list(seq_len(num_leaves)), as.character(initial_tree$frame$yval[leaves]))

      # Extract design matrix for evalrand
      preds_evalrand <- predict(initial_tree, rdf_evalrand)
      nodes_of_evalrand <- unlist(lookup_preds[as.character(preds_evalrand)])
      nis_evalrand <- table(rdf_evalrand$id)
      I.evalrand <- length(nis_evalrand)
      epsilon_evalrand <- rdf_evalrand$y - preds_evalrand

      # Optimise rho (or keep fixed if fixrho a fixed constant)
      if (!fixrho) {
        rho_optim <- optimise_randeff(num_leaves, I.evalrand, nis_evalrand, nodes_of_evalrand, epsilon_evalrand,
                                    x0dataset=x0dataset, tree=initial_tree, lookup_preds=lookup_preds,
                                    test_data=test_data)
      } else {
        if (!is.numeric(fixrho) || length(fixrho) != 1 || fixrho < -1 || fixrho > 1) stop("Error: Parameter must be a numeric constant between -1 and 1")
        rho_optim <- fixrho
      }
      # Extract design matrix for evalfix
      preds_evalfix <- predict(initial_tree, rdf_evalfix)
      nodes_of_evalfix <- unlist(lookup_preds[as.character(preds_evalfix)])
      nis_evalfix <- table(rdf_evalfix$id)
      I.evalfix <- length(nis_evalfix)

      XWX_XWY_mats <- XWX_XWY_calc(rho_optim, num_leaves, I.evalfix, nis_evalfix, nodes_of_evalfix, rdf_evalfix$y)
      final_tree <- initial_tree
      final_tree$frame$yval[final_tree$frame$var=="<leaf>"] <- solve(XWX_XWY_mats[[1]]+diag(1e-6,num_leaves,num_leaves), XWX_XWY_mats[[2]])
      forest[[b]] <- final_tree
    }
    forests[[l]] <- forest
  }
  endtime <- proc.time()
  if (floor((endtime - starttime)[[3]]/60)!=0) print(paste0("Total Runtime: ", floor((endtime - starttime)[[3]]/60)," minutes" ))

  num_variables <- length(attr(terms(formula), "term.labels"))
  res <- list(forest=forests, call=match.call(), num.samples=nrow(data), num.groups=length(unique(data$id)), num.of.covariates=num_variables, beta=beta, honesty=honesty, weight_optimiser=weight_optimiser, correlation=correlation)
  class(res$forest) <- "crf.forest"
  class(res) <- "crf"
  return(res)
}
