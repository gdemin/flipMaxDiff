#' \code{LatentClassMaxDiff}
#' @description Fits a latent class rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
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
#' @export
LatentClassMaxDiff <- function(design, version, best, worst, alternative.names, n.classes,
                               subset = NULL, weights = NULL, seed = 123)
{
    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights)
    respondent.indices <- dat$respondent.indices
    n.respondents <- length(respondent.indices)
    n.beta <- dat$n.alternatives - 1
    n.tasks <- dat$n.tasks
    X <- dat$X
    weights <- dat$weights
    tol <- 0.0000001

    class.memberships <- randomClassMemberships(n.respondents, n.classes, seed)
    p <- inferParameters(class.memberships, X, weights, respondent.indices, n.beta)

    # These are the
    # p <- c(1.361,
    #        0.5023,
    #        1.8976,
    #        1.1576,
    #        2.001,
    #        1.7042,
    #        0.8686,
    #        0.3861,
    #        0.8134,
    #        -2.386,
    #        -3.5968,
    #        -1.8937,
    #        -3.1371,
    #        -2.3911,
    #        -2.5816,
    #        -3.4425,
    #        -3.6203,
    #        -3.5637,
    #        -0.3314001)
    previous.ll <- -Inf

    repeat
    {
        pp <- posteriorProbabilities(p, X, weights, respondent.indices, n.classes, n.beta)
        p <- inferParameters(pp, X, weights, respondent.indices, n.beta)
        ll <- logLikelihood(p, X, weights, respondent.indices, n.classes, n.beta)
        # print(ll)
        if (ll - previous.ll < tol)
            break
        else
            previous.ll <- ll
    }
    list(log.likelihood = ll, coef = p)
}

# Randomly assign n individuals to classes
randomClassMemberships <- function(n, n.classes, seed = 123)
{
    set.seed(seed)
    memberships <- vector("integer", n)
    repeat
    {
        memberships <- sample(1:n.classes, n, replace = T)
        # Ensure that each class has at least one member
        if (length(table(memberships)) == n.classes)
            break
    }
    res <- matrix(0, n, n.classes)
    for (i in 1:n)
        res[i, memberships[i]] <- 1
    res
}

# Infer parameters given class memberships
inferParameters <- function(class.memberships, X, weights, respondent.indices, n.beta)
{
    n.classes <- ncol(class.memberships)
    n.parameters <- n.beta * n.classes + n.classes - 1
    res <- vector("numeric", n.parameters)
    # Class parameters
    for (c in 1:n.classes)
    {
        class.weights <- vector("numeric", nrow(X))
        for (l in 1:length(respondent.indices))
            class.weights[respondent.indices[[l]]] <- class.memberships[l, c]
        class.weights <- class.weights * weights
        solution <- optimizeMaxDiff(X, class.weights, n.beta)
        res[((c - 1) * n.beta + 1):(c * n.beta)] <- solution$par
    }

    # Class size parameters
    sizes <- colMeans(class.memberships)
    res[(n.classes * n.beta + 1):n.parameters] <- -log(sizes[2:n.classes] / sizes[1])
    res
}

# Calculate posterior probabilities of class membership given parameters
posteriorProbabilities <- function(pars, X, weights, respondent.indices, n.classes, n.beta)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, weights, respondent.indices, n.classes)

    n.repsondents <- length(respondent.indices)
    repsondents.log.densities <- vector("numeric", n.repsondents)
    for (l in 1:n.repsondents)
        repsondents.log.densities[l] <- logOfSum(log.class.weights + log.densities[l, ])


    exp(t(matrix(rep(log.class.weights, n.repsondents), n.classes)) + log.densities
        - t(matrix(rep(repsondents.log.densities, each = n.classes), n.classes)))
}

logLikelihood <- function(pars, X, weights, respondent.indices, n.classes, n.beta)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, weights, respondent.indices, n.classes)
    res <- 0
    for (l in 1:length(respondent.indices))
        res <- res + logOfSum(log.class.weights + log.densities[l, ])
    res
}

logDensities <- function(pars.list, X, weights, respondent.indices, n.classes)
{
    log.shares <- matrix(NA, nrow(X), n.classes)
    for (c in 1:n.classes)
    {
        e.u <- matrix(exp(c(0, pars.list[[c]])[X]), ncol = ncol(X))
        log.shares[, c] <- logDensitiesBestWorst(e.u, weights)
    }

    res <- matrix(NA, length(respondent.indices), n.classes)
    for (l in 1:length(respondent.indices))
        res[l, ] <- colSums(log.shares[respondent.indices[[l]], , drop = FALSE])
    res
}

# A numerically stable way of calculating log(sum(exp(logs)))
logOfSum <- function(logs)
{
    maxlog <- max(logs)
    maxlog + log(sum(exp(logs - maxlog)))
}

# Extract class sizes from parameters
getClassWeights <- function(pars, n.classes, n.beta)
{
    exp.pars <- exp(-pars[(n.classes * n.beta + 1):length(pars)])
    total <- 1 + sum(exp.pars)
    c(1, exp.pars) / total
}

# Extract parameters for each class
getClassParameters <- function(pars, n.classes, n.beta)
{
    pars.list <- vector("list", n.classes)
    for (c in 1:n.classes)
        pars.list[[c]] <- pars[((c - 1) * n.beta + 1):(c * n.beta)]
    pars.list
}
