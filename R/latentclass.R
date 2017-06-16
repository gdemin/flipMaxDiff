latentClassMaxDiff <- function(dat, ind.levels, resp.pars = NULL, n.classes = 1, seed = 123,
                               initial.parameters = NULL, n.previous.parameters = 0, trace = 0,
                               apply.weights = TRUE, tol = 0.0001)
{
    n.respondents <- length(dat$respondent.indices)
    n.levels <- length(ind.levels)
    n.beta <- dat$n.alternatives - 1
    n.questions <- dat$n.questions.in
    n.choices <- ncol(dat$X.in)
    X <- dat$X.in
    weights <- dat$weights
    alternative.names <- dat$alternative.names

    boost <- computeBoost(resp.pars, X, n.respondents, n.questions, n.choices)

    p <- if (!is.null(initial.parameters))
        initial.parameters
    else
    {
        class.memberships <- randomClassMemberships(n.levels, n.classes, seed)
        inferParameters(class.memberships, X, boost, weights, ind.levels, n.beta, trace)
    }

    previous.ll <- -Inf

    repeat
    {
        pp <- posteriorProbabilities(p, X, boost, ind.levels, n.classes, n.beta)
        p <- inferParameters(pp, X, boost, weights, ind.levels, n.beta, trace)
        ll <- logLikelihood(p, X, boost, weights, ind.levels, n.classes, n.beta, n.questions, apply.weights = apply.weights)
        # print(ll)
        if (ll - previous.ll < tol)
            break
        else
            previous.ll <- ll
    }

    respondent.pp <- matrix(NA, n.respondents * n.questions, n.classes)
    for (l in 1:n.levels)
    {
        n.ind <- length(ind.levels[[l]])
        respondent.pp[ind.levels[[l]], ] <- t(matrix(rep(pp[l, ], n.ind), ncol = n.ind))
    }
    respondent.pp <- respondent.pp[(1:n.respondents) * n.questions, , drop = FALSE]
    if (!is.null(dat$subset))
    {
        posterior.probabilities <- matrix(NA, length(dat$subset), n.classes)
        posterior.probabilities[dat$subset, ] <- respondent.pp
    }
    else
        posterior.probabilities <- respondent.pp

    if (!is.null(dat$subset) && !is.null(resp.pars))
    {
        input.respondent.pars <- matrix(NA, length(dat$subset), ncol(resp.pars))
        input.respondent.pars[dat$subset, ] <- resp.pars
    }
    else
        input.respondent.pars <- resp.pars

    result <- list(posterior.probabilities = posterior.probabilities,
                   log.likelihood = ll,
                   n.classes = n.classes)
    if (n.classes > 1)
    {
        coef <- matrix(0, nrow = n.beta + 1, ncol = n.classes)
        for (c in 1:n.classes)
            coef[2:(n.beta + 1), c] <- p[((c - 1) * n.beta + 1):(c * n.beta)]
        rownames(coef) <- dat$alternative.names
        colnames(coef) <- paste("Class", 1:n.classes)
        result$class.sizes <- getClassWeights(p, n.classes, n.beta)
        result$class.preference.shares <- classPreferenceShares(coef, result$class.sizes)
    }
    else
    {
        coef <- c(0, p)
        names(coef) <- dat$alternative.names
        result$class.sizes <- 1
        result$class.preference.shares <- exp(coef) / sum(exp(coef))
    }
    result$coef <- coef
    result$effective.sample.size <- ess <- sum(weights) / n.questions
    result$n.parameters <- numberOfParameters(n.beta, n.classes, "Finite", FALSE) + n.previous.parameters
    result$bic <- -2 * ll + log(ess) * result$n.parameters
    result$respondent.parameters <- computeRespondentParameters(result, alternative.names, input.respondent.pars)
    class(result) <- "FitMaxDiff"
    result
}

# Randomly assign n individuals to classes
randomClassMemberships <- function(n, n.classes, seed = 123)
{
    set.seed(seed)
    memberships <- vector("integer", n)
    if (n < n.classes)
        stop("The number of individuals may not be less than the number of classes.")

    # Could be optimized to avoid the repeat loop which is slow when n = n.classes and both are large.
    # Do this by assigning one individual to a class at the start.
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
inferParameters <- function(class.memberships, X, boost, weights, ind.levels, n.beta, trace)
{
    n.classes <- ncol(class.memberships)
    n.parameters <- n.beta * n.classes + n.classes - 1
    n.levels <- length(ind.levels)
    res <- vector("numeric", n.parameters)
    # Class parameters
    for (c in 1:n.classes)
    {
        class.weights <- vector("numeric", nrow(X))
        for (l in 1:n.levels)
            class.weights[ind.levels[[l]]] <- class.memberships[l, c]
        if (!is.null(weights))
            class.weights <- class.weights * weights
        solution <- optimizeMaxDiff(X, boost, class.weights, n.beta, trace)
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
posteriorProbabilities <- function(pars, X, boost, ind.levels, n.classes, n.beta)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, boost, ind.levels, n.classes)

    n.levels <- length(ind.levels)
    repsondents.log.densities <- vector("numeric", n.levels)
    for (l in 1:n.levels)
        repsondents.log.densities[l] <- logOfSum(log.class.weights + log.densities[l, ])

    exp(t(matrix(rep(log.class.weights, n.levels), n.classes)) + log.densities
        - t(matrix(rep(repsondents.log.densities, each = n.classes), n.classes)))
}

logLikelihood <- function(pars, X, boost, weights, ind.levels, n.classes, n.beta, n.questions, apply.weights)
{
    log.class.weights <- log(getClassWeights(pars, n.classes, n.beta))
    class.pars <- getClassParameters(pars, n.classes, n.beta)
    log.densities <- logDensities(class.pars, X, boost, ind.levels, n.classes)
    n.levels <- length(ind.levels)
    res <- 0
    for (l in 1:n.levels)
    {
        if (apply.weights)
        {
            res <- res + logOfSum(log.class.weights + log.densities[l, ]) * weights[(l - 1) * n.questions + 1]
        }
        else
        {
            res <- res + logOfSum(log.class.weights + log.densities[l, ])
        }
    }
    res
}

logDensities <- function(pars.list, X, boost, ind.levels, n.classes)
{
    log.shares <- matrix(NA, nrow(X), n.classes)
    for (c in 1:n.classes)
    {
        e.u <- exp(matrix(c(0, pars.list[[c]])[X], ncol = ncol(X)) + boost)
        log.shares[, c] <- logDensitiesBestWorst(e.u, rep(1, length(e.u)))
    }

    n.levels <- length(ind.levels)
    res <- matrix(NA, n.levels, n.classes)
    for (l in 1:n.levels)
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

getLevelIndices <- function(characteristic, n.questions)
{
    n.respondents <- length(characteristic)
    lvls <- unique(characteristic)
    n.levels <- length(lvls)
    result <- vector("list", n.levels)
    for (l in 1:n.levels)
        result[[l]] <- (1:(n.respondents * n.questions))[rep(characteristic == lvls[l], each = n.questions)]
    result
}

classPreferenceShares <- function(coef, class.sizes)
{
    n.classes <- ncol(coef)
    pref.shares <- exp(coef) / t(matrix(rep(colSums(exp(coef)), nrow(coef)), nrow = n.classes))
    result <- matrix(NA, nrow(pref.shares), ncol(pref.shares) + 1)
    colnames(result) <- c(paste("Class", 1:n.classes), "Total")
    rownames(result) <- rownames(pref.shares)
    result[, 1:ncol(pref.shares)] <- pref.shares
    result[, ncol(pref.shares) + 1] <- rowSums(pref.shares *
        t(matrix(rep(class.sizes, nrow(pref.shares)), ncol = nrow(pref.shares))))
    result
}

computeRespondentParameters <- function(object, alternative.names, input.respondent.pars = NULL)
{
    pp <- object$posterior.probabilities
    coef <- object$coef
    n.classes <- object$n.classes
    if (n.classes > 1)
        result <- pp[ , 1, drop = FALSE] %*% t(coef[, 1, drop = FALSE])
    else
        result <- pp[ , 1, drop = FALSE] %*% t(coef)
    if (n.classes > 1)
        for (c in 2:n.classes)
            result <- result + pp[ , c, drop = FALSE] %*% t(coef[, c, drop = FALSE])
    if (!is.null(input.respondent.pars))
        result <- result + input.respondent.pars

    colnames(result) <- alternative.names
    result
}

computeBoost <- function(resp.pars, X, n.respondents, n.questions, n.choices)
{
    if (is.null(resp.pars))
        result <- matrix(0, nrow = n.respondents * n.questions, ncol = n.choices)
    else
    {
        resp.pars <- as.matrix(resp.pars)
        result <- matrix(0, nrow = n.respondents * n.questions, ncol = n.choices)
        for (i in 1:n.respondents)
        {
            ind <- ((i - 1) * n.questions + 1):(i * n.questions)
            result[ind, ] <- matrix(resp.pars[i, X[ind, ]], ncol = n.choices)
        }
    }
    result
}
