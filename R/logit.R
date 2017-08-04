#' @importFrom stats cor optim
optimizeMaxDiff <- function(X, boost, weights, n.pars, trace = 0, is.tricked = FALSE)
{
    init.b <- seq(.01, .02, length.out = n.pars)
    optim(init.b,
          logLikelihoodMaxDiff,
          gr = gradientMaxDiffR,
          # gr = NULL,
          X = X,
          boost = boost,
          weights = weights,
          is.tricked = is.tricked,
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
#' @param is.tricked Whether to use tricked logit instead of rank-ordered logit with ties.
logLikelihoodMaxDiff <- function(b, X, boost, weights, is.tricked)
{
    b[b > 100] = 100
    b[b < -100] = -100

    e.u <- exp(matrix(c(0, b)[X], ncol = ncol(X)) + boost)
    logDensityMaxDiff(e.u, weights, is.tricked)
}

gradientMaxDiffR <- function(b, X, boost, weights, is.tricked)
{
    b[b > 100] = 100
    b[b < -100] = -100
    e.u <- exp(matrix(c(0, b)[X], ncol = ncol(X)) + boost)
    gradientMaxDiff(e.u, X - 1, weights, length(b), is.tricked)
}
