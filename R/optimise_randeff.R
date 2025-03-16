# Optimise for random effects
#' @importFrom stats optim

optimise_randeff_ptwvar_equicorr <- function(num_leaves, I, nis, nodesis, epsilon, x0dataset, tree, lookup_preds, ...) {
  predict_x0 <- predict(tree, x0dataset)
  eRx <- rep(0,num_leaves)
  eRx[lookup_preds[as.character(predict_x0)][[1]]] <- 1
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_equicorr_cpp(rho, num_leaves, I, nis, nodesis, epsilon)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- t(eRx) %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX %*% eRx
    return(loss)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
  return(rho_optim)
}

optimise_randeff_ptwvar_ar1 <- function(num_leaves, I, nis, nodesis, epsilon, x0dataset, tree, lookup_preds, ...) {
  predict_x0 <- predict(tree, x0dataset)
  eRx <- rep(0,num_leaves)
  eRx[lookup_preds[as.character(predict_x0)][[1]]] <- 1
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_ar1_cpp(rho, num_leaves, I, nis, nodesis, epsilon)
    solveXWX_eRx <- solve( XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves) , eRx)
    loss <- t(solveXWX_eRx) %*% XWX_XWSWX_mats[[2]] %*% solveXWX_eRx
    return(loss)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
  return(rho_optim)
}

optimise_randeff_trainmse_equicorr <- function(num_leaves, I, nis, nodesis, epsilon, ...) {
  # x0, tree, convert_pred_to_leaf defunct
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XX_equicorr_cpp(rho, num_leaves, I, nis, nodesis, epsilon)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
}

optimise_randeff_trainmse_ar1 <- function(num_leaves, I, nis, nodesis, epsilon, ...) {
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XX_ar1_cpp(rho, num_leaves, I, nis, nodesis, epsilon)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
}
optimise_randeff_testmse_equicorr <- function(num_leaves, I, nis, nodesis, epsilon, test_data, test_density, tree, lookup_preds, ...) {
  if (is.null(test_data) && !is.null(test_density)) stop("test_density not currently supported. Create your own dataset following the chosen test_density and input it in test_data")
  N.TEST <- nrow(test_data)
  preds_test <- predict(tree, test_data)
  nodesisTEST <- unlist(lookup_preds[as.character(preds_test)])
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XtestX_equicorr_cpp(rho, num_leaves, I, nis, nodesis, epsilon, N.TEST, nodesisTEST)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss/N.TEST)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
  return(rho_optim)
}

optimise_randeff_testmse_ar1 <- function(num_leaves, I, nis, nodesis, epsilon, test_data, test_density, tree, lookup_preds, ...) {
  if (is.null(test_data) && !is.null(test_density)) stop("test_density not currently supported. Create your own dataset following the chosen test_density and input it in test_data")
  N.TEST <- nrow(test_data)
  preds_test <- predict(tree, test_data)
  nodesisTEST <- unlist(lookup_preds[as.character(preds_test)])
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XtestX_ar1_cpp(rho, num_leaves, I, nis, nodesis, epsilon, N.TEST, nodesisTEST)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss/N.TEST)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
  return(rho_optim)
}
