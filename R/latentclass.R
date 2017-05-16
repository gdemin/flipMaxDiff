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
inferParameters <- function(class.memberships, X, weights, ind.levels, n.beta, trace)
{
    n.classes <- ncol(class.memberships)
    n.parameters <- n.beta * n.classes + n.classes - 1
    res <- vector("numeric", n.parameters)
    # Class parameters
    for (c in 1:n.classes)
    {
        class.weights <- vector("numeric", nrow(X))
        for (l in 1:length(ind.levels))
            class.weights[ind.levels[[l]]] <- class.memberships[l, c]
        class.weights <- class.weights * weights
        solution <- optimizeMaxDiff(X, class.weights, n.beta, trace)
        res[((c - 1) * n.beta + 1):(c * n.beta)] <- solution$par
    }

    # Class size parameters
    if (n.classes > 1)
    {
        sizes <- colMeans(class.memberships)
        res[(n.classes * n.beta + 1):n.parameters] <- -log(sizes[2:n.classes] / sizes[1])
    }
    res
}

# Calculate posterior probabilities of class membership given parameters
posteriorProbabilities <- function(pars, X, ind.levels, n.classes, n.beta)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, ind.levels, n.classes)

    n.repsondents <- length(ind.levels)
    repsondents.log.densities <- vector("numeric", n.repsondents)
    for (l in 1:n.repsondents)
        repsondents.log.densities[l] <- logOfSum(log.class.weights + log.densities[l, ])


    exp(t(matrix(rep(log.class.weights, n.repsondents), n.classes)) + log.densities
        - t(matrix(rep(repsondents.log.densities, each = n.classes), n.classes)))
}

logLikelihood <- function(pars, X, weights, ind.levels, n.classes, n.beta)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, ind.levels, n.classes)
    n.tasks <- length(weights) / length(ind.levels)
    res <- 0
    for (l in 1:length(ind.levels))
        res <- res + logOfSum(log.class.weights + log.densities[l, ]) * weights[(l - 1) * n.tasks + 1]
    res
}

logDensities <- function(pars.list, X, ind.levels, n.classes)
{
    log.shares <- matrix(NA, nrow(X), n.classes)
    for (c in 1:n.classes)
    {
        e.u <- matrix(exp(c(0, pars.list[[c]])[X]), ncol = ncol(X))
        log.shares[, c] <- logDensitiesBestWorst(e.u, rep(1, length(e.u)))
    }

    res <- matrix(NA, length(ind.levels), n.classes)
    for (l in 1:length(ind.levels))
        res[l, ] <- colSums(log.shares[ind.levels[[l]], , drop = FALSE])
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
    if (n.classes > 1)
    {
        exp.pars <- exp(-pars[(n.classes * n.beta + 1):length(pars)])
        total <- 1 + sum(exp.pars)
        c(1, exp.pars) / total
    }
    else
        1
}

# Extract parameters for each class
getClassParameters <- function(pars, n.classes, n.beta)
{
    pars.list <- vector("list", n.classes)
    for (c in 1:n.classes)
        pars.list[[c]] <- pars[((c - 1) * n.beta + 1):(c * n.beta)]
    pars.list
}
