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
#' @param sub.model.outputs If TRUE, prints diagnostics on interrum models.
#' @param trace Non-negative integer indicating the detail of outputs provided when fitting models: 0 indicates
#' @param lc Whether to run latent class step at the end if characteristics are supplied.
#' no outputs, and 6 is the most detailed outputs.
#' @param output Output type. Can be "Probabilities" or "Classes".
#' @export
FitMaxDiff <- function(design, version, best, worst, alternative.names, n.classes = 1,
                       subset = NULL, weights = NULL, characteristics = NULL, seed = 123,
                       initial.parameters = NULL, trace = 0, sub.model.outputs = FALSE, lc = TRUE, output = "Probabilities")
{
    if (!is.null(weights) && !is.null(characteristics))
        stop("Weights are not able to be applied when characteristics are supplied")
    if (!lc && is.null(characteristics))
        stop("There is no model to run. Select covariates and/or run latent class analysis over respondents.")

    apply.weights <- is.null(characteristics)

    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights, characteristics)
    n.respondents <- length(dat$respondent.indices)
    n.tasks <- dat$n.tasks
    resp.pars <- NULL
    n.previous.parameters <- 0
    characteristics <- dat$characteristics
    covariates.notes <- NULL
    if (!is.null(characteristics))
    {
        resp.pars <- matrix(0, nrow = n.respondents, ncol = dat$n.alternatives)
        n.characteristics <- length(characteristics)
        covariates.notes <- character(n.characteristics)
        for (i in 1:n.characteristics)
        {
            ind.levels <- getLevelIndices(characteristics[[i]], n.tasks)
            n.levels <- length(ind.levels)
            best.bic <- Inf
            best.solution <- NULL
            # Loop over all possible class sizes
            for (n.c in 1:n.levels)
            {
                solution <- latentClassMaxDiff(dat, alternative.names, ind.levels, resp.pars, n.c, seed,
                                   initial.parameters, n.previous.parameters, trace, apply.weights = apply.weights)
                if (solution$bic < best.bic)
                {
                    best.bic <- solution$bic
                    best.solution <- solution
                }
            }
            if (best.solution$n.classes > 1)
            {
                resp.pars <- as.matrix(RespondentParameters(best.solution))
                n.previous.parameters <- best.solution$n.parameters
                covariates.notes[i] <- paste0(names(characteristics)[i], " - ", best.solution$n.classes, " classes")
            }
            else
            {
                covariates.notes[i] <- paste0(names(characteristics)[i], " - no selection")
                if (sub.model.outputs)
                    cat("Covariate:", names(characteristics)[i],
                        "BIC:", best.solution$bic,
                        "LL:", best.solution$log.likelihood,
                        "Number of classes:", best.solution$n.classes,
                        (if (best.solution$n.classes == 1) "Excluded" else "Included in final model"),
                        "\n")
            }
        }
    }

    result <- if (lc || is.null(characteristics))
        latentClassMaxDiff(dat, alternative.names, dat$respondent.indices, resp.pars, n.classes, seed,
                       initial.parameters, n.previous.parameters, trace, apply.weights = apply.weights)
    else
        best.solution
    if (sub.model.outputs)
        cat("Latent class analysis",
            "BIC:", result$bic,
            "LL:", result$log.likelihood,
            "Number of classes:", result$n.classes)
    result$n.respondents <- n.respondents
    result$subset <- subset
    result$weights <- weights
    result$n.tasks <- n.tasks
    result$n.alternatives.per.task <- ncol(dat$X)
    result$covariates.notes <- covariates.notes
    result$output <- output
    result$lc <- lc
    result
}

latentClassMaxDiff <- function(dat, alternative.names, ind.levels, resp.pars = NULL, n.classes = 1, seed = 123,
                               initial.parameters = NULL, n.previous.parameters = 0, trace = 0, apply.weights = TRUE)
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
        inferParameters(class.memberships, X, boost, weights, ind.levels, n.beta, trace)
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
        result$class.sizes <- getClassWeights(p, n.classes, n.beta)
        probabilities <- exp(coef) / t(matrix(rep(colSums(exp(coef)), nrow(coef)), nrow = n.classes))
        table.data <- matrix(NA, nrow(probabilities), ncol(probabilities) + 1)
        colnames(table.data) <- c(paste("Class", 1:n.classes), "Total")
        rownames(table.data) <- rownames(probabilities)
        table.data[, 1:ncol(probabilities)] <- probabilities
        table.data[, ncol(probabilities) + 1] <- rowSums(probabilities * t(matrix(rep(result$class.sizes, nrow(probabilities)), ncol = nrow(probabilities))))
        probabilities <- table.data
    }
    else
    {
        coef <- c(0, p)
        n.parameters <- length(p)
        names(coef) <- dat$alternative.names
        probabilities <- exp(coef) / sum(exp(coef))
        result$class.sizes <- 1
    }
    result$coef <- coef
    result$effective.sample.size <- ess <- sum(weights) / n.tasks
    result$n.parameters <- n.parameters + n.previous.parameters
    result$bic <- -2*ll + log(ess) * result$n.parameters
    result$probabilities <- probabilities
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
        result <- pp[ , 1, drop = FALSE] %*% t(coef[, 1, drop = FALSE])
    else
        result <- pp[ , 1, drop = FALSE] %*% t(coef)
    if (n.classes > 1)
        for (c in 2:n.classes)
            result <- result + pp[ , c, drop = FALSE] %*% t(coef[, c, drop = FALSE])
    if (!is.null(object$input.respondent.pars))
        result <- result + object$input.respondent.pars
    result <- as.data.frame(result)

    names(result) <- if (n.classes > 1) rownames(coef) else names(coef)
    result
}

#' \code{Memberships}
#' @description Allocates respondents to classes based on the maximum posterior probabilty.
#' @param object A \code{FitMaxDiff} object.
#' @export
Memberships <- function(object)
{
    pp <- object$posterior.probabilities
    .fun <- function(x) match(max(x), x)[1]
    apply(pp, 1, .fun)
}

#' \code{print.CorrespondenceAnalysis}
#' @param x FitMaxDiff object.
#' @param ... further arguments passed to or from other methods.
#' @importFrom flipFormat MaxDiffTable MaxDiffTableClasses FormatAsPercent FormatWithDecimals
#' @importFrom stats median quantile
#' @export
print.FitMaxDiff <- function(x, ...)
{
    title <- "Max-Diff: Latent Class Analysis"
    footer <- paste0("n = ", x$n.respondents, "; ")
    if (!is.null(x$subset))
        footer <- paste0(footer, "Filters have been applied; ")
    if (!is.null(x$weights))
        footer <- paste0(footer, "Weights have been applied; Effective sample size: ",
                         FormatWithDecimals(x$effective.sample.size, 2), "; ")
    footer <- paste0(footer, "Number of questions: ", x$n.tasks, "; ")
    footer <- paste0(footer, "Alternatives per question: ", x$n.alternatives.per.task, "; ")
    footer <- paste0(footer, "Log-likelihood: ", FormatWithDecimals(x$log.likelihood, 2), "; ")
    footer <- paste0(footer, "BIC: ", FormatWithDecimals(x$bic, 2), "; ")
    if (!x$lc)
        footer <- paste0(footer, "Latent class analysis over respondents not applied; ")

    if (x$n.classes == 1 && is.null(x$covariates.notes))
    {
        col.labels <- "Probabilities (%)"
        MaxDiffTableClasses(as.matrix(x$probabilities), col.labels, title, "", footer)
    }
    else if (x$output == "Probabilities")
    {
        subtitle <- if (!is.null(x$covariates.notes))
            paste0("Covariates: ", paste(x$covariates.notes, collapse = ", "))
        else
            "Covariates: none"

        # subtitle with covariates
        resp.pars <- as.matrix(RespondentParameters(x))
        probs <- exp(resp.pars) / rowSums(exp(resp.pars))
        stats.table <- matrix(NA, nrow = ncol(probs), ncol = 6)
        for (i in 1:ncol(probs))
        {
            p <- probs[, i]
            stats.table[i, 1] <- mean(p)
            stats.table[i, 2] <- median(p)
            stats.table[i, 3] <- quantile(p, 0.25)
            stats.table[i, 4] <- quantile(p, 0.75)
            stats.table[i, 5] <- min(p)
            stats.table[i, 6] <- max(p)
        }
        rownames(stats.table) <- colnames(probs)
        MaxDiffTable(stats.table, title, subtitle, footer)
    }
    else if (x$output == "Classes")
    {
        if (!is.null(x$covariates.notes))
            stop("Class table cannot be displayed when covariates are applied.")
        col.labels <- c(paste("Class", 1:x$n.classes, "(%)<br>Size:", FormatAsPercent(x$class.sizes, 3)), "Total")
        MaxDiffTableClasses(as.matrix(x$probabilities), col.labels, title, "", footer)
    }
}
