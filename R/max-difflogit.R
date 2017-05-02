#' \code{dMaxDiff}
#' @description The density of a respondent's max-diff data.
#' @param b A vector of parameter estimates.
#' @param x The experimental design for a respondent  (a \code{\link{list}}).
dMaxDiff <- function(b, x)
{
    prob <- prod(as.numeric(lapply(x, function(x.i) dBestWorst(exp(b[x.i])))))
    prob

}


#' \code{dMaxDiff}
#' @description The log-likelihood for a max-diff experiment.
#' @param b A vector of parameter estimates.
#' @param X The experimental design for a sample (a \code{\link{list}})
logLikelihoodMaxDiff = function(b, X)
{
   b[b > 100] = 100
   b[b < -100] = -100
   probs <- as.numeric(lapply(X, b = c(0, b), dMaxDiff))
   sum(log(probs))
}


#' \code{FitMaxDiff}
#' @description Fits a rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst A matrix of integers showing the choice of 'worst'.
#' @param names The names of the alternatives.
#' @export
FitMaxDiff <- function(design, version, best, worst, names)
{
    X <- IntegrateDesignAndData(design, version, best, worst)
    n.alternatives <- max(design[, -1:-2])
    init.b <- seq(.01,.02, length.out = n.alternatives - 1)
    solution = optim(init.b, logLikelihoodMaxDiff,  gr = NULL, X = X, method =  "BFGS", control = list(fnscale  = -1, maxit = 1000, trace = 6), hessian = FALSE)
    pars = c(0, solution$par)
    names(pars) = names
    #names(pars) = dimnames(stacked.data)[[2]]
    list(log.likelihood = solution$value, coef = pars)
}




