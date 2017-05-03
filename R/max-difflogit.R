#' \code{FitMaxDiff}
#' @description Fits a rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst A matrix of integers showing the choice of 'worst'.
#' @param names A character vector names of the alternatives. If only a single element is supplied, it is split by commas.
#' @param subset An optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights An optional vector of sampling or frequency weights.
#' @param trace Non-negative integer indicating the detail of outputs provided during estimation: 0 indicates
#' no outputs, and 6 is the most detailed outputs.
#' @importFrom flipData CalibrateWeight CleanSubset CleanWeights
#' @export
FitMaxDiff <- function(design, version, best, worst, names, subset = NULL, weights = NULL, trace = 0)
{
    # Cleaning and checking data
    n <- length(best[[1]])
    if (!is.null(weights))
        weights <- CleanWeights(weights)
    subset <- CleanSubset(subset, n)
    if (is.null(version))
        version <- rep(1, n)
    version <- version[subset]
    if (!is.null(weights))
    {
        weights <- weights[subset]
        weights <- CalibrateWeight(weights)
    }
    best <- as.data.frame(best[subset, ])
    worst <- as.data.frame(worst[subset, ])
    if (!("Version" %in% names(design)))
        design <- cbind(Version = 1, design)
    if (!("Task" %in% names(design)))
        design <- cbind(Task = 1:nrow(design), design)
    X <- IntegrateDesignAndData(design, version, best, worst)
    n.alternatives <- max(design[, -1:-2])
    if (length(names) != n.alternatives)
        names <- strsplit(names)
    if (length(names) != n.alternatives)
        stop("The number of 'names' does not match the number of alternatives in the 'design'.")
    # Estimating
    init.b <- seq(.01,.02, length.out = n.alternatives - 1)
    solution = optim(init.b, logLikelihoodMaxDiff,
                     gr = NULL,
                     X = X,
                     weights = weights,
                     method =  "BFGS",
                     control = list(fnscale  = -1, maxit = 1000, trace = trace), hessian = FALSE)
    pars = c(0, solution$par)
    names(pars) = names
    list(log.likelihood = solution$value, coef = pars, n = n)
}


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
logLikelihoodMaxDiff = function(b, X, weights)
{
   b[b > 100] = 100
   b[b < -100] = -100
   probs <- as.numeric(lapply(X, b = c(0, b), dMaxDiff))
   log.probs <- log(probs)
   if (!is.null(weights))
       log.probs <- log.probs * weights
   sum(log.probs)
}


