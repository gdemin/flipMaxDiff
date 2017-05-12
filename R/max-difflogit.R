#' \code{FitMaxDiff}
#' @description Fits a rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst A matrix of integers showing the choice of 'worst'.
#' @param alternative.names A character vector names of the alternatives. If only a single element is supplied, it is split by commas.
#' @param subset An optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights An optional vector of sampling or frequency weights.
#' @param trace Non-negative integer indicating the detail of outputs provided during estimation: 0 indicates
#' no outputs, and 6 is the most detailed outputs.
#' @importFrom flipData CalibrateWeight CleanSubset CleanWeights
#' @importFrom flipFormat TrimWhitespace
#' @importFrom stats cor optim
#' @export
FitMaxDiff <- function(design, version, best, worst, alternative.names, subset = NULL, weights = NULL, trace = 0)
{
    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights)
    solution <- optimizeMaxDiff(dat$X, dat$weights, dat$n.alternatives - 1, trace)
    pars = c(0, solution$par)
    names(pars) = dat$alternative.names
    list(log.likelihood = solution$value, coef = pars, n = dat$n)
}

cleanAndCheckData <- function(design, version, best, worst, alternative.names, subset = NULL, weights = NULL)
{
    # Cleaning and checking data
    n <- length(best[[1]])
    n.tasks <- nrow(design)
    if (missing(version))
        version <- rep(1, n)
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
    weights <- if (is.null(weights))
        rep(1, n.tasks * n)
    else
        rep(weights, each = n.tasks)
    best <- as.data.frame(best[subset, ])
    worst <- as.data.frame(worst[subset, ])
    if (!("Version" %in% names(design)))
        design <- cbind(Version = 1, design)
    versions.in.variable <- sort(unique(version))
    versions.in.design <- sort(unique(design$Version))
    compare.version <- all.equal(versions.in.variable, versions.in.design)
    if (is.character(compare.version))
        stop("The 'design' and 'version' have incompatible version numbers.")
    if (!("Task" %in% names(design)))
        design <- cbind(Task = 1:n.tasks, design)
    dat <- IntegrateDesignAndData(design, version, best, worst)
    n.alternatives <- max(design[, -1:-2])
    if (missing(alternative.names))
        alternative.names <- paste("Alternative", 1:n.alternatives)
    if (length(alternative.names) != n.alternatives)
    {
        alternative.names <- strsplit(alternative.names, ",")
        alternative.names <- TrimWhitespace(alternative.names)
    }
    if (length(alternative.names) != n.alternatives)
        stop("The number of 'alternative.names' does not match the number of alternatives in the 'design'.")

    list(X = dat$X,
         weights = weights,
         alternative.names = alternative.names,
         n = n,
         n.alternatives = n.alternatives,
         n.tasks = n.tasks,
         respondent.indices = dat$respondent.indices)
}

optimizeMaxDiff <- function(X, weights, n.pars, trace = 0)
{
    init.b <- seq(.01, .02, length.out = n.pars)
    optim(init.b,
          logLikelihoodMaxDiff,
          gr = gradientMaxDiff,
          X = X,
          weights = weights,
          method =  "BFGS",
          control = list(fnscale  = -1, maxit = 1000, trace = trace),
          hessian = FALSE)
}

#' \code{dMaxDiff}
#' @description The log-likelihood for a max-diff experiment.
#' @param b A vector of parameter estimates.
#' @param X The experimental design for a sample (a \code{\link{list}})
#' @param weights An optional vector of sampling or frequency weights.
logLikelihoodMaxDiff = function(b, X, weights)
{
    b[b > 100] = 100
    b[b < -100] = -100
    e.u <- matrix(exp(c(0, b)[X]), ncol = ncol(X))
    logDensityBestWorst(e.u, weights)
}

gradientMaxDiff = function(b, X, weights)
{
    b[b > 100] = 100
    b[b < -100] = -100
    e.u <- matrix(exp(c(0, b)[X]), ncol = ncol(X))
    gradientBestWorst(e.u, X - 1, weights, length(b))
}
