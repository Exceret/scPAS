####################################################
#####  Augmented and Penalized Minimization    #####
#####  APML0 (L1/L2/Laplacian+L0)              #####
#####  Penalty: L0, L1, L2, Laplacian          #####
#####  Algorithm: one-step coordinate descent  #####
####################################################

#' @title Augmented and Penalized Minimization with L0 Regularization
#' @description
#' Fits regularized regression models with L0 penalty for feature selection.
#' Supports Gaussian, binomial, and Cox regression with various penalty types
#' including Lasso, Elastic Net, and Network regularization.
#'
#' @param x Input matrix of features (n x p)
#' @param y Response vector/matrix:
#' - For Gaussian: numeric vector
#' - For binomial: binary vector (0/1)
#' - For Cox: matrix with columns "time" and "status"
#' @param family Type of regression model: "gaussian", "binomial", or "cox"
#' @param penalty Type of penalty: "Lasso" (L1), "Enet" (L1 + L2), or "Net" (L1 + Laplacian)
#' @param Omega Optional penalty matrix for network regularization (for "Net" penalty)
#' @param alpha Elastic net mixing parameter (0 <= alpha <= 1).
#' alpha=1 for Lasso, alpha=0 for Ridge
#' @param lambda Regularization parameter sequence. If NULL, automatically generated
#' @param nlambda Number of lambda values to generate (default: 50)
#' @param rlambda Ratio of smallest to largest lambda (default: adaptive)
#' @param wbeta Adaptive weights for coefficients (default: all 1)
#' @param sgn Sign constraints for coefficients (default: all 1)
#' @param nfolds Number of folds for cross-validation (default: 1)
#' @param foldid Optional vector specifying fold membership for each observation
#' @param ill Whether to use log-likelihood for model selection (default: TRUE)
#' @param iL0 Whether to use L0 regularization (default: TRUE)
#' @param icutB Whether to use coefficient thresholding (default: FALSE)
#' @param ncutB Number of coefficient cuts for thresholding (default: 10)
#' @param ifast Whether to use fast algorithm (default: TRUE)
#' @param isd Whether to standardize features (default: FALSE)
#' @param iysd Whether to standardize response (Gaussian only, default: FALSE)
#' @param ifastr Whether to use fast cross-validation (default: TRUE)
#' @param keep.beta Whether to keep all coefficient paths (default: FALSE)
#' @param thresh Convergence threshold for optimization (default: 1e-6)
#' @param maxit Maximum number of iterations (default: 1e+5)
#' @param threshC Convergence threshold for coordinate descent (default: 1e-5)
#' @param maxitC Maximum iterations for coordinate descent (default: 1e+2)
#' @param threshP Convergence threshold for proximal operator (default: 1e-5)
#'
#' @return An object of class "APML0" containing:
#' \item{Beta}{Estimated coefficients}
#' \item{fit}{Model fit information including lambda path and performance metrics}
#' \item{family}{Type of regression model}
#' \item{penalty}{Type of penalty used}
#' \item{adaptive}{Whether adaptive penalties were used}
#' \item{flag}{Convergence flags}
#'
#'
#' @examples
#' \dontrun{
#' # Gaussian regression example
#' set.seed(123)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#'
#' # Fit Lasso model with L0 regularization
#' fit <- APML0(x, y, family = "gaussian", penalty = "Lasso", nlambda = 30)
#'
#' # Cox regression example
#' library(survival)
#' y_surv <- cbind(time = rexp(100), status = rbinom(100, 1, 0.5))
#' fit_cox <- APML0(x, y_surv, family = "cox", penalty = "Enet", alpha = 0.5)
#' }
#'
#' @seealso
#' [print.APML0()] for printing results, [LmL0()], [LogL0()], [CoxL0()]
#' @export
APML0 = function(
  x,
  y,
  family = c("gaussian", "binomial", "cox"),
  penalty = c("Lasso", "Enet", "Net"),
  Omega = NULL,
  alpha = 1.0,
  lambda = NULL,
  nlambda = 50,
  rlambda = NULL,
  wbeta = rep(1, ncol(x)),
  sgn = rep(1, ncol(x)),
  nfolds = 1,
  foldid = NULL,
  ill = TRUE,
  iL0 = TRUE,
  icutB = FALSE,
  ncutB = 10,
  ifast = TRUE,
  isd = FALSE,
  iysd = FALSE,
  ifastr = TRUE,
  keep.beta = FALSE,
  thresh = 1e-6,
  maxit = 1e+5,
  threshC = 1e-5,
  maxitC = 1e+2,
  threshP = 1e-5
) {
  #fcall=match.call()
  family = match.arg(family)
  penalty = match.arg(penalty)

  if (penalty == "Net" & is.null(Omega)) {
    penalty = "Enet"
    cat("Enet was performed as no input of Omega")
  }
  if (penalty %in% c("Enet", "Net") & alpha == 1.0) {
    penalty = "Lasso"
    cat("Lasso was performed as alpha=1.0")
  }

  if (alpha != 1.0) {
    if (is.null(Omega)) {
      penalty = "Enet"
    } else if (!is.null(Omega)) {
      penalty = "Net"
    }
  } else {
    penalty = "Lasso"
  }

  wbeta = abs(wbeta)

  fit = switch(
    family,
    "gaussian" = LmL0(
      x,
      y,
      Omega,
      alpha,
      lambda,
      nlambda,
      rlambda,
      wbeta,
      sgn,
      nfolds,
      foldid,
      ill,
      iL0,
      icutB,
      ncutB,
      ifast,
      isd,
      iysd,
      keep.beta,
      thresh,
      maxit
    ),
    "binomial" = LogL0(
      x,
      y,
      Omega,
      alpha,
      lambda,
      nlambda,
      rlambda,
      wbeta,
      sgn,
      nfolds,
      foldid,
      iL0,
      icutB,
      ncutB,
      ifast,
      isd,
      keep.beta,
      thresh,
      maxit,
      threshC,
      maxitC,
      threshP
    ),
    "cox" = CoxL0(
      x,
      y,
      Omega,
      alpha,
      lambda,
      nlambda,
      rlambda,
      wbeta,
      sgn,
      nfolds,
      foldid,
      iL0,
      icutB,
      ncutB,
      ifast,
      isd,
      ifastr,
      keep.beta,
      thresh,
      maxit
    )
  )
  fit$family = family

  #fit$call=fcall
  class(fit) = "APML0"
  return(fit)
}
