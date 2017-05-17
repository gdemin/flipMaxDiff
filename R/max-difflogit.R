#' \code{FitMaxDiff}
#' @description Fits a latent class rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each task.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst A matrix of integers showing the choice of 'worst'.
#' @param alternative.names A character vector names of the alternatives. If only a single element is supplied, it is split by commas.
#' @param n.classes The number of latent classes.
#' @param subset An optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights An optional vector of sampling or frequency weights.
#' @param characteristics Data frame of characteristics on which to run varying coefficients by latent class boosting.
#' @param seed Seed for initial random class assignments.
#' @param initial.parameters Specify initial parameters intead of starting at random.
#' @param trace Non-negative integer indicating the detail of outputs provided during estimation: 0 indicates
#' @param lc Whether to run latent class step at the end if characteristics are supplied.
#' no outputs, and 6 is the most detailed outputs.
#' @export
FitMaxDiff <- function(design, version, best, worst, alternative.names, n.classes = 1,
                       subset = NULL, weights = NULL, characteristics = NULL, seed = 123,
                       initial.parameters = NULL, trace = 0, lc = TRUE)
{
    if (!is.null(weights) && !is.null(characteristics))
        stop("Weights are not able to be applied when characteristics are supplied")

    apply.weights <- is.null(characteristics)

    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights, characteristics)
    n.respondents <- length(dat$respondent.indices)
    n.tasks <- dat$n.tasks
    resp.pars <- NULL
    characteristics <- dat$characteristics
    if (!is.null(characteristics))
    {
        resp.pars <- matrix(0, nrow = n.respondents, ncol = dat$n.alternatives)
        for (ch in characteristics)
        {
            ind.levels <- getLevelIndices(ch, n.tasks)
            n.levels <- length(ind.levels)
            best.bic <- Inf
            best.solution <- NULL
            # Loop over all possible class sizes
            for (n.c in 1:n.levels)
            {
                solution <- latentClassMaxDiff(dat, alternative.names, ind.levels, resp.pars, n.c, seed,
                                   initial.parameters, trace, apply.weights = apply.weights)
                if (solution$bic < best.bic)
                {
                    best.bic <- solution$bic
                    best.solution <- solution
                }
            }
            resp.pars <- as.matrix(RespondentParameters(best.solution))
        }
    }
    if (lc || is.null(characteristics))
        latentClassMaxDiff(dat, alternative.names, dat$respondent.indices, resp.pars, n.classes, seed,
                       initial.parameters, trace, apply.weights = apply.weights)
    else
        best.solution
}

latentClassMaxDiff <- function(dat, alternative.names, ind.levels, resp.pars = NULL, n.classes = 1, seed = 123,
                               initial.parameters = NULL, trace = 0, apply.weights = TRUE)
{
    n.respondents <- length(dat$respondent.indices)
    n.levels <- length(ind.levels)
    n.beta <- dat$n.alternatives - 1
    n.tasks <- dat$n.tasks
    n.choices <- ncol(dat$X)
    X <- dat$X
    weights <- dat$weights
    tol <- 0.00001

    boost <- if (is.null(resp.pars))
        matrix(0, nrow = n.respondents * n.tasks, ncol = n.choices)
    else
    {
        resp.pars <- as.matrix(resp.pars)
        tmp <- matrix(0, nrow = n.respondents * n.tasks, ncol = n.choices)
        for (i in 1:n.respondents)
        {
            ind <- ((i - 1) * n.tasks + 1):(i * n.tasks)
            tmp[ind, ] <- matrix(resp.pars[i, X[ind, ]], ncol = n.choices)
        }
        tmp
    }


    p <- if (!is.null(initial.parameters))
        initial.parameters
    else
    {
        class.memberships <- randomClassMemberships(n.levels, n.classes, seed)
        p <- inferParameters(class.memberships, X, boost, weights, ind.levels, n.beta, trace)
    }

    previous.ll <- -Inf

    repeat
    {
        pp <- posteriorProbabilities(p, X, boost, ind.levels, n.classes, n.beta)
        p <- inferParameters(pp, X, boost, weights, ind.levels, n.beta, trace)
        ll <- logLikelihood(p, X, boost, weights, ind.levels, n.classes, n.beta, n.tasks, apply.weights = apply.weights)
        # print(ll)
        if (ll - previous.ll < tol)
            break
        else
            previous.ll <- ll
    }

    respondent.pp <- matrix(NA, n.respondents * n.tasks, n.classes)
    for (l in 1:n.levels)
    {
        n.ind <- length(ind.levels[[l]])
        respondent.pp[ind.levels[[l]], ] <- t(matrix(rep(pp[l, ], n.ind), ncol = n.ind))
    }
    respondent.pp <- respondent.pp[(1:n.respondents) * n.tasks, , drop = FALSE]

    result <- list(posterior.probabilities = respondent.pp,
                   log.likelihood = ll,
                   n.classes = n.classes,
                   input.respondent.pars = resp.pars)
    if (n.classes > 1)
    {
        coef <- matrix(0, nrow = n.beta + 1, ncol = n.classes)
        for (c in 1:n.classes)
            coef[2:(n.beta + 1), c] <- p[((c - 1) * n.beta + 1):(c * n.beta)]
        rownames(coef) <- dat$alternative.names
        colnames(coef) <- paste("Class", 1:n.classes)
        n.parameters <- prod(dim(coef[-1,])) + n.classes - 1
        result$class.sizes = getClassWeights(p, n.classes, n.beta)
    }
    else
    {
        coef <- c(0, p)
        n.parameters <- length(p)
        names(coef) <- dat$alternative.names
    }
    result$coef = coef
    result$effective.sample.size <- ess <- sum(weights) / n.tasks
    result$bic = -2*ll + log(ess) * n.parameters
    class(result) <- "FitMaxDiff"
    result
}

#' @importFrom flipData CalibrateWeight CleanSubset CleanWeights
#' @importFrom flipFormat TrimWhitespace
cleanAndCheckData <- function(design, version, best, worst, alternative.names, subset = NULL, weights = NULL,
                              characteristics = NULL)
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
        rep(1, n.tasks * length(version))
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
        alternative.names <- strsplit(alternative.names, ",")[[1]]
        alternative.names <- TrimWhitespace(alternative.names)
    }
    if (length(alternative.names) != n.alternatives)
        stop("The number of 'alternative.names' does not match the number of alternatives in the 'design'.")
    if (!is.null(characteristics))
        characteristics <- characteristics[subset, , drop = FALSE]

    list(X = dat$X,
         weights = weights,
         alternative.names = alternative.names,
         n = n,
         n.alternatives = n.alternatives,
         n.tasks = n.tasks,
         respondent.indices = dat$respondent.indices,
         characteristics = characteristics)
}

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

#' \code{RespondentParameters}
#' @description Computes parameters for each respondent.
#' @param object A \code{FitMaxDiff} object.
#' @export
RespondentParameters <- function(object)
{
    pp <- object$posterior.probabilities
    coef <- object$coef
    n.classes <- object$n.classes
    if (n.classes > 1)
        result <- pp[ ,1, drop = FALSE] %*% t(coef[, 1, drop = FALSE])
    else
        result <- pp[ ,1, drop = FALSE] %*% t(coef)
    if (n.classes > 1)
        for (c in 2:n.classes)
            result <- result + pp[ ,c , drop = FALSE] %*% t(coef[, c, drop = FALSE])
    result <- result + object$input.respondent.pars
    result <- as.data.frame(result)
    names(result) <- rownames(coef)
    result
}
