#' @importFrom MASS mvrnorm
mixtureOfNormalsMaxDiff <- function(dat, n.classes, normal.covariance, seed = 123, initial.parameters = NULL,
                                    trace = 0, pool.variance = FALSE, n.draws = 100, is.tricked = FALSE)
{
    n.respondents <- length(dat$respondent.indices)
    n.beta <- dat$n.alternatives - 1
    n.questions <- dat$n.questions.in
    max.iterations.since.best <- 20
    X <- dat$X.in
    weights <- dat$weights
    alternative.names <- dat$alternative.names

    if (is.null(initial.parameters))
    {
        class.memberships <- randomClassMemberships(n.respondents, n.classes, seed)
        class.weights <- colMeans(class.memberships)
        means <- generateInitialMeans(class.memberships, X, weights, n.questions, n.beta, trace, is.tricked)
        covariances <- generateInitialCovariances(n.classes, n.beta)
    }
    else
    {
        class.weights <- initial.parameters$class.weights
        means <- initial.parameters$means
        covariances <- initial.parameters$covariances
    }

    best.log.likelihood <- -Inf
    best.pars <- list(class.weights = class.weights, means = means, covariances = covariances)
    iterations.since.best <- 0
    repeat
    {
        log.likelihood <- 0
        sum.weights <- vector("numeric", n.classes)
        sum.weighted.beta <- vector("list", n.classes)
        sum.weighted.squares <- vector("list", n.classes)
        beta.class.draws <- vector("list", n.classes)
        for (c in 1:n.classes)
        {
            sum.weighted.beta[[c]] <- vector("numeric", n.beta + 1)
            sum.weighted.squares[[c]] <- matrix(0, n.beta + 1, n.beta + 1)
            beta.class.draws[[c]] <- mvrnorm(n.draws * n.respondents, means[[c]], covariances[[c]])
        }

        # For each respondent, calculate the log choice probability
        for (i in 1:n.respondents)
        {
            beta.draws <- matrix(NA, n.classes * n.draws, n.beta + 1)
            for (c in 1:n.classes)
                beta.draws[(c - 1) * n.draws + (1:n.draws), ] <- beta.class.draws[[c]][(i - 1) * n.draws + (1:n.draws), ]
            ind <- (i - 1) * n.questions + (1:n.questions)

            bd <- beta.draws[, 2:(n.beta + 1)] - beta.draws[, 1]

            logs.of.k <- logKernels(bd, X[ind, ], rep(1, n.questions), is.tricked)
            log.choice.p <- logChoiceProb(logs.of.k, class.weights)

            resp.weight <- weights[ind[1]]
            log.likelihood <- log.likelihood + log.choice.p * resp.weight

            draw.weights <- exp(rep(log(class.weights), each = n.draws) + logs.of.k - log.choice.p)
            for (c in 1:n.classes)
            {
                ind.draws <- (c - 1) * n.draws + (1:n.draws)
                w <- draw.weights[ind.draws]
                b <- beta.draws[ind.draws, ]
                sum.weights[c] <- sum.weights[c] + sum(w) * resp.weight
                sum.weighted.beta[[c]] <- sum.weighted.beta[[c]] + colSums(b * w) * resp.weight
                sum.weighted.squares[[c]] <- sum.weighted.squares[[c]] + sumWeightedOuterProducts(b, w) * resp.weight
            }
        }

        class.weights <- prop.table(sum.weights)
        means <- parameterMeans(sum.weights, sum.weighted.beta)
        covariances <- parameterCovariances(means, sum.weights, sum.weighted.beta, sum.weighted.squares,
                                            pool.variance, normal.covariance)

        # print(log.likelihood)

        if (log.likelihood > best.log.likelihood)
        {
            best.log.likelihood <- log.likelihood
            best.pars <- list(class.weights = class.weights, means = means, covariances = covariances)
            iterations.since.best <- 0
        } else
            iterations.since.best <- iterations.since.best + 1

        if (iterations.since.best == max.iterations.since.best)
            break
    }

    info <- respondentLevelInfo(X, best.pars$class.weights, best.pars$means, best.pars$covariances,
                                alternative.names, dat$subset, n.classes, n.respondents, n.questions,
                                1000, n.beta, is.tricked)

    best.covariances <- lapply(best.pars$covariances, function(x) {
        # rownames(x) <- alternative.names[-1]
        # colnames(x) <- alternative.names[-1]
        rownames(x) <- alternative.names
        colnames(x) <- alternative.names
        x
    })

    result <- list(log.likelihood = best.log.likelihood)
    if (n.classes > 1)
    {
        coef <- matrix(0, nrow = n.beta + 1, ncol = n.classes)
        for (c in 1:n.classes)
            coef[, c] <- best.pars$means[[c]]
        rownames(coef) <- dat$alternative.names
        colnames(coef) <- paste("Class", 1:n.classes)
        result$class.sizes <- best.pars$class.weights
        result$class.preference.shares <- classPreferenceShares(coef, best.pars$class.weights)
    }
    else
    {
        # coef <- c(0, best.pars$means[[1]])
        coef <- best.pars$means[[1]]
        names(coef) <- dat$alternative.names
        result$class.sizes <- 1
        result$class.preference.shares <- exp(coef) / sum(exp(coef))
    }
    result$coef <- coef
    result$covariances <- best.covariances
    result$effective.sample.size <- ess <- sum(weights) / n.questions
    result$n.parameters <- numberOfParameters(n.beta, n.classes, TRUE, normal.covariance, pool.variance)
    result$bic <- -2 * best.log.likelihood + log(ess) * result$n.parameters
    result$posterior.probabilities <- info$posterior.probabilities
    result$respondent.parameters <- info$respondent.parameters
    class(result) <- "FitMaxDiff"
    result
}

generateInitialMeans <- function(class.memberships, X, weights, n.questions, n.beta, trace, is.tricked)
{
    n.classes <- ncol(class.memberships)
    n.respondents <- nrow(class.memberships)
    res <- vector("list", n.classes)
    boost <- rep(0, nrow(X))
    for(c in 1:n.classes)
    {
        w <- weights * rep(class.memberships[, c], each = n.questions)
        res[[c]] <- optimizeMaxDiff(X, boost, w, n.beta, trace, is.tricked)$par
        res[[c]] <- c(0, res[[c]])
    }
    res
}

generateInitialCovariances <- function(n.classes, n.beta)
{
    res <- vector("list", n.classes)
    for(c in 1:n.classes)
        res[[c]] <- diag(n.beta + 1) * 0.1
    res
}

parameterMeans <- function(sum.weights, sum.weighted.beta)
{
    n.classes <- length(sum.weights)
    res <- vector("list", n.classes)
    for (c in 1:n.classes)
    {
        res[[c]] <- sum.weighted.beta[[c]] / sum.weights[c]
        res[[c]][1] <- 0
    }
    res
}

parameterCovariances <- function(means, sum.weights, sum.weighted.beta, sum.weighted.squares, pool.variance,
                                 normal.covariance)
{
    n.classes <- length(sum.weights)
    n.beta <- length(sum.weighted.beta[[1]])
    res <- vector("list", n.classes)
    if (pool.variance)
    {
        pooled.covariance <- matrix(0, n.beta, n.beta)
        for (c in 1:n.classes)
        {
            cross.term <- outer(means[[c]], sum.weighted.beta[[c]])
            pooled.covariance <- pooled.covariance + sum.weighted.squares[[c]] - cross.term - t(cross.term) +
                outer(means[[c]], means[[c]]) * sum.weights[c]
        }
        pooled.covariance <- constrainCovariance(pooled.covariance / sum(sum.weights), normal.covariance)
        for (c in 1:n.classes)
            res[[c]] <- pooled.covariance
    }
    else
    {
        for (c in 1:n.classes)
        {
            cross.term <- outer(means[[c]], sum.weighted.beta[[c]])
            covariance.matrix <- (sum.weighted.squares[[c]] - cross.term - t(cross.term) +
                                      outer(means[[c]], means[[c]]) * sum.weights[c]) / sum.weights[c]
            res[[c]] <- constrainCovariance(covariance.matrix, normal.covariance)
        }
    }
    res
}

logChoiceProb <- function(logs.of.k, class.weights)
{
    n.classes <- length(class.weights)
    n.draws <- length(logs.of.k) / n.classes
    logs.class.sum <- vector("numeric", n.classes)
    for (c in 1:n.classes)
        logs.class.sum[c] <- log(class.weights[c]) + logOfSum(logs.of.k[(c - 1) * n.draws + (1:n.draws)]) - log(n.draws)
    logOfSum(logs.class.sum)
}

numberOfParameters <- function(n.beta, n.classes, is.mixture.of.normals = FALSE, normal.covariance = NULL,
                               pool.variance = FALSE)
{
    result <- n.beta * n.classes + n.classes - 1
    if (is.mixture.of.normals)
    {
        if (normal.covariance == "Spherical")
        {
            if (pool.variance)
                result <- result + 1
            else
                result <- result + n.classes
        }
        else if (normal.covariance == "Diagonal")
        {
            if (pool.variance)
                result <- result + n.beta
            else
                result <- result + n.beta * n.classes
        }
        else if (normal.covariance == "Full")
        {
            if (pool.variance)
                result <- result + 0.5 * n.beta * (n.beta + 1)
            else
                result <- result + 0.5 * n.beta * (n.beta + 1) * n.classes
        }
        else
            stop(paste("Distribution not handled:", normal.covariance))
    }
    result
}

constrainCovariance <- function(covariance.matrix, normal.covariance)
{
    if (normal.covariance == "Spherical")
        result <- diag(rep(mean(diag(covariance.matrix)), ncol(covariance.matrix)))
    else if (normal.covariance == "Diagonal")
        result <- diag(diag(covariance.matrix))
    else if (normal.covariance == "Full")
        result <- covariance.matrix
    else
        stop(paste("Distribution not handled:", normal.covariance))
    result
}

respondentLevelInfo <- function(X, class.weights, means, covariances, alternative.names, subset, n.classes,
                                n.respondents, n.questions, n.draws, n.beta, is.tricked)
{
    beta.class.draws <- vector("list", n.classes)
    for (c in 1:n.classes)
        beta.class.draws[[c]] <- mvrnorm(n.draws * n.respondents, means[[c]], covariances[[c]])

    pp <- matrix(NA, n.respondents, n.classes)
    resp.pars <- matrix(0, n.respondents, n.beta + 1)
    for (i in 1:n.respondents)
    {
        beta.draws <- matrix(NA, n.classes * n.draws, n.beta + 1)
        for (c in 1:n.classes)
            beta.draws[(c - 1) * n.draws + (1:n.draws), ] <- beta.class.draws[[c]][(i - 1) * n.draws + (1:n.draws), ]
        ind <- (i - 1) * n.questions + (1:n.questions)

        bd <- beta.draws[, 2:(n.beta + 1)] - beta.draws[, 1]

        logs.of.k <- logKernels(bd, X[ind, ], rep(1, n.questions), is.tricked)
        log.choice.p <- logChoiceProb(logs.of.k, class.weights)
        draw.weights <- exp(rep(log(class.weights), each = n.draws) + logs.of.k - log.choice.p)

        pp[i, ] <- prop.table(colSums(matrix(draw.weights, nrow = n.draws)))
        # resp.pars[i, 2:(n.beta + 1)] <- colSums(beta.draws * draw.weights) / sum(draw.weights)
        resp.pars[i, ] <- colSums(beta.draws * draw.weights) / sum(draw.weights)
    }
    posterior.probabilities <- matrix(NA, length(subset), n.classes)
    posterior.probabilities[subset, ] <- pp
    respondent.parameters <- matrix(NA, length(subset), n.beta + 1)
    respondent.parameters[subset, ] <- resp.pars

    colnames(respondent.parameters) <- alternative.names
    list(posterior.probabilities = posterior.probabilities, respondent.parameters = respondent.parameters)
}
