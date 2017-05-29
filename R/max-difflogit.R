#' \code{FitMaxDiff}
#' @description Fits a latent class rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task' or 'Question',
#' and the remaining variables contain the alternatives shown in each task.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A data frame of factors or a matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to correspond to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst As with 'best', except denoting worst..
#' @param worst A matrix of integers showing the choice of 'worst'.
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
#' @param tasks.left.out Number of tasks to leave out for cross-validation.
#' @export
FitMaxDiff <- function(design, version = NULL, best, worst, alternative.names, n.classes = 1,
                       subset = NULL, weights = NULL, characteristics = NULL, seed = 123,
                       initial.parameters = NULL, trace = 0, sub.model.outputs = FALSE, lc = TRUE,
                       output = "Probabilities", tasks.left.out = 0)
{
    if (!is.null(weights) && !is.null(characteristics))
        stop("Weights are not able to be applied when characteristics are supplied")
    if (!lc && is.null(characteristics))
        stop("There is no model to run. Select covariates and/or run latent class analysis over respondents.")

    apply.weights <- is.null(characteristics)

    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights,
                             characteristics, seed, tasks.left.out)
    n.respondents <- length(dat$respondent.indices)
    n.questions.in <- dat$n.questions.in
    resp.pars <- NULL
    n.previous.parameters <- 0
    characteristics <- dat$characteristics
    covariates.notes <- NULL
    covariates.chosen <- FALSE
    if (!is.null(characteristics))
    {
        resp.pars <- matrix(0, nrow = n.respondents, ncol = dat$n.alternatives)
        n.characteristics <- length(characteristics)
        covariates.notes <- character(n.characteristics)
        for (i in 1:n.characteristics)
        {
            ind.levels <- getLevelIndices(characteristics[[i]], n.questions.in)
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
                resp.pars <- as.matrix(RespondentParameters(best.solution))[dat$subset, ]
                n.previous.parameters <- best.solution$n.parameters
                covariates.notes[i] <- paste0(names(characteristics)[i], " - ", best.solution$n.classes, " classes")
                covariates.chosen <- TRUE
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
    else if (covariates.chosen)
        best.solution
    else
        stop("No model applied. Choose different covariates or enable latent class analysis over respondents.")

    result$in.sample.accuracy <- predictionAccuracy(result, dat$X.in, n.questions.in, dat$subset)
    result$out.sample.accuracy <- if (tasks.left.out > 0)
        predictionAccuracy(result, dat$X.out, tasks.left.out, dat$subset)
    else
        NA

    if (sub.model.outputs)
        cat("Latent class analysis",
            "BIC:", result$bic,
            "LL:", result$log.likelihood,
            "Number of classes:", result$n.classes)
    result$n.respondents <- n.respondents # this should be the unfiltered number of respondents
    result$subset <- subset
    result$weights <- weights
    result$n.questions <- dat$n.questions
    result$n.alternatives.per.task <- ncol(dat$X.in)
    result$covariates.notes <- covariates.notes
    result$output <- output
    result$lc <- lc
    result$tasks.left.out <- tasks.left.out

    resp.pars <- as.matrix(RespondentParameters(result))[dat$subset, ]
    result$respondent.probabilities <- exp(resp.pars) / rowSums(exp(resp.pars))

    result
}

latentClassMaxDiff <- function(dat, alternative.names, ind.levels, resp.pars = NULL, n.classes = 1, seed = 123,
                               initial.parameters = NULL, n.previous.parameters = 0, trace = 0, apply.weights = TRUE)
{
    n.respondents <- length(dat$respondent.indices)
    n.levels <- length(ind.levels)
    n.beta <- dat$n.alternatives - 1
    n.questions <- dat$n.questions.in
    n.choices <- ncol(dat$X.in)
    X <- dat$X.in
    weights <- dat$weights
    tol <- 0.0001
    boost <- if (is.null(resp.pars))
        matrix(0, nrow = n.respondents * n.questions, ncol = n.choices)
    else
    {
        resp.pars <- as.matrix(resp.pars)
        tmp <- matrix(0, nrow = n.respondents * n.questions, ncol = n.choices)
        for (i in 1:n.respondents)
        {
            ind <- ((i - 1) * n.questions + 1):(i * n.questions)
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
                   n.classes = n.classes,
                   input.respondent.pars = input.respondent.pars)
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
    result$effective.sample.size <- ess <- sum(weights) / n.questions
    result$n.parameters <- n.parameters + n.previous.parameters
    result$bic <- -2*ll + log(ess) * result$n.parameters
    result$probabilities <- probabilities
    class(result) <- "FitMaxDiff"
    result
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

predictionAccuracy <- function(object, X, n.questions, subset)
{
    score <- rep(NA, nrow(X))
    resp.pars <- as.matrix(RespondentParameters(object)[subset, ])
    for (i in 1:nrow(resp.pars))
    {
        pars <- resp.pars[i, ]
        for (j in 1:n.questions)
        {
            ind <- (i - 1) * n.questions + j
            u <- pars[X[ind, ]]
            score[ind] <- all(u[1] > u[-1])
        }
    }
    sum(score) / nrow(X)
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
    .fun <- function(x)
    {
        if (any(is.na(x)))
            NA
        else
            match(max(x), x)[1]
    }
    apply(pp, 1, .fun)
}

#' \code{print.CorrespondenceAnalysis}
#' @param x FitMaxDiff object.
#' @param ... further arguments passed to or from other methods.
#' @importFrom flipFormat HistTable MaxDiffTableClasses FormatAsPercent FormatWithDecimals
#' @importFrom stats median quantile
#' @export
print.FitMaxDiff <- function(x, ...)
{
    title <- if (!is.null(x$covariates.notes))
        "Max-Diff: Varying Coefficients"
    else
        "Max-Diff: Latent Class Analysis"
    footer <- paste0("n = ", x$n.respondents, "; ")
    if (!is.null(x$subset) && !all(x$subset))
        footer <- paste0(footer, "Filters have been applied; ")
    if (!is.null(x$weights))
        footer <- paste0(footer, "Weights have been applied; Effective sample size: ",
                         FormatWithDecimals(x$effective.sample.size, 2), "; ")
    footer <- paste0(footer, "Number of questions: ", x$n.questions, "; ")
    if (x$tasks.left.out > 0)
    {
        footer <- paste0(footer, "Questions used in estimation: ", x$n.questions - x$tasks.left.out, "; ")
        footer <- paste0(footer, "Questions left out: ", x$tasks.left.out, "; ")
    }
    footer <- paste0(footer, "Alternatives per question: ", x$n.alternatives.per.task, "; ")
    footer <- paste0(footer, "Log-likelihood: ", FormatWithDecimals(x$log.likelihood, 2), "; ")
    footer <- paste0(footer, "BIC: ", FormatWithDecimals(x$bic, 2), "; ")
    footer <- if (!x$lc)
        paste0(footer, "Latent class analysis over respondents not applied; ")
    else
    {
        if (x$n.classes == 1)
            paste0(footer, "Latent class analysis over respondents: ", x$n.classes, " class; ")
        else
            paste0(footer, "Latent class analysis over respondents: ", x$n.classes, " classes; ")
    }

    subtitle <- if (!is.na(x$out.sample.accuracy))
        paste0("Prediction accuracy (leave-", x$tasks.left.out , "-out cross-validation): ",
               FormatAsPercent(x$out.sample.accuracy, 3))
    else
        paste0("Prediction accuracy (in-sample): ", FormatAsPercent(x$in.sample.accuracy, 3))

    if (x$n.classes == 1 && is.null(x$covariates.notes))
    {
        col.labels <- "Probabilities (%)"
        MaxDiffTableClasses(as.matrix(x$probabilities), col.labels, title, subtitle, footer)
    }
    else if (x$output == "Probabilities")
    {
        if (!is.null(x$covariates.notes))
            subtitle <- c(subtitle, paste0("Covariates: ", paste(x$covariates.notes, collapse = ", ")))

        probs <- x$respondent.probabilities
        stats.table <- matrix(NA, nrow = ncol(probs), ncol = 1)
        for (i in 1:ncol(probs))
            stats.table[i, 1] <- FormatWithDecimals(mean(probs[, i], na.rm = TRUE) * 100, 1)
        colnames(stats.table) <- "Mean Probability (%)"

        HistTable(100 * probs, title = title, subtitle = subtitle, footer = footer,
                  bin.size = 5, bin.min = 0, bin.max = 100, hist.width = 100, hist.height = 20, stats.table)
    }
    else if (x$output == "Classes")
    {
        if (!is.null(x$covariates.notes))
            stop("Class table cannot be displayed when covariates are applied.")
        col.labels <- c(paste("Class", 1:x$n.classes, "(%)<br>Size:", FormatAsPercent(x$class.sizes, 3)), "Total")
        MaxDiffTableClasses(as.matrix(x$probabilities), col.labels, title, subtitle, footer)
    }
}
