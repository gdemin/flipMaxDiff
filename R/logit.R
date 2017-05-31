#' @importFrom stats cor optim
optimizeMaxDiff <- function(X, boost, weights, n.pars, trace = 0)
{
    init.b <- seq(.01, .02, length.out = n.pars)
    optim(init.b,
          logLikelihoodMaxDiff,
          gr = gradientMaxDiff,
          X = X,
          boost = boost,
          weights = weights,
          method =  "BFGS",
          control = list(fnscale  = -1, maxit = 1000, trace = trace),
          hessian = FALSE)
}

#' \code{dMaxDiff}
#' @description The log-likelihood for a max-diff experiment.
#' @param b A vector of parameter estimates.
#' @param X The experimental design for a sample (a \code{\link{list}})
#' @param boost Boosting constants.
#' @param weights An optional vector of sampling or frequency weights.
logLikelihoodMaxDiff <- function(b, X, boost, weights)
{
    b[b > 100] = 100
    b[b < -100] = -100

    e.u <- exp(matrix(c(0, b)[X], ncol = ncol(X)) + boost)
    logDensityBestWorst(e.u, weights)
}

gradientMaxDiff <- function(b, X, boost, weights)
{
    b[b > 100] = 100
    b[b < -100] = -100
    e.u <- exp(matrix(c(0, b)[X], ncol = ncol(X)) + boost)
    gradientBestWorst(e.u, X - 1, weights, length(b))
}
