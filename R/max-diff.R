#' \code{FitMaxDiff}
#' @description Fits a latent class rank-ordered logit model with ties to a max-diff experiment.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task' or 'Question',
#' and the remaining variables contain the alternatives shown in each task.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A data frame of factors or a matrix of integers showing the choices made by each respondent on each of the questions. One column
#' for each task. The integers need to correspond to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst As with 'best', except denoting worst..
#' @param alternative.names A \code{character} vector of the alternative names. Where \code{best}
#' and  \code{worst} are factors or characters, these names must match them.
#' @param n.classes The number of latent classes.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights An optional vector of sampling or frequency weights.
#' @param characteristics Data frame of characteristics on which to run varying coefficients by latent class boosting.
#' @param seed Seed for initial random class assignments.
#' @param initial.parameters Specify initial parameters intead of starting at random.
#' @param sub.model.outputs If TRUE, prints diagnostics on interim models.
#' @param trace Non-negative integer indicating the detail of outputs provided when fitting models: 0 indicates
#' no outputs, and 6 is the most detailed outputs.
#' @param lc Whether to run latent class step at the end if characteristics are supplied.
#' @param output Output type. Can be "Probabilities" or "Classes".
#' @param tasks.left.out Number of questions to leave out for cross-validation.
#' @param algorithm If "HB", Hierarchical Bayes with a MVN prior is used and the other parameters are ignored.
#' @param is.mixture.of.normals Whether to model with mixture of normals instead of LCA.
#' @param normal.covariance The form of the covariance matrix for mixture of normals.
#' Can be 'Full, 'Spherical', 'Diagonal'.
#' @param pool.variance Whether to pool parameter covariances between classes in mixture of normals.
#' @param lc.tolerance The tolerance used for defining convergence in latent class analysis.
#' @param n.draws The number of draws when fitting mixture of normals.
#' @param is.tricked Whether to use tricked logit instead of rank-ordered logit with ties.
#' @param hb.iterations The number of iterations in Hierarchical Bayes.
#' @param hb.chains The number of chains in Hierarchical Bayes.
#' @param hb.max.tree.depth http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#' @param hb.adapt.delta http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#' @export
FitMaxDiff <- function(design, version = NULL, best, worst, alternative.names, n.classes = 1,
                       subset = NULL, weights = NULL, characteristics = NULL, seed = 123,
                       initial.parameters = NULL, trace = 0, sub.model.outputs = FALSE, lc = TRUE,
                       output = "Probabilities", tasks.left.out = 0, is.mixture.of.normals = FALSE,
                       algorithm = "Default", normal.covariance = "Full", pool.variance = FALSE,
                       lc.tolerance = 0.0001, n.draws = 100, is.tricked = FALSE,
                       hb.iterations = 100, hb.chains = 1, hb.max.tree.depth = 10,
                       hb.adapt.delta = 0.8)
{
    if (!is.null(weights) && !is.null(characteristics))
        stop("Weights are not able to be applied when characteristics are supplied.")
    if (!is.null(weights) && algorithm == "HB")
        stop("Weights are not able to be applied for Hierarchical Bayes.")
    if (!lc && is.null(characteristics))
        stop("There is no model to run. Select covariates and/or run latent class analysis over respondents.")
    if (!is.null(characteristics) && is.mixture.of.normals)
        stop("Mixture of normals cannot be selected when characteristics are supplied.")

    apply.weights <- is.null(characteristics)
    questions.left.out <- tasks.left.out # we now refer to tasks as questions

    dat <- cleanAndCheckData(design, version, best, worst, alternative.names, subset, weights,
                             characteristics, seed, questions.left.out)

    if (algorithm == "HB")
    {
        result <- hierarchicalBayesMaxDiff(dat, hb.iterations, hb.chains, hb.max.tree.depth,
                                           hb.adapt.delta, TRUE)
    }
    else if (is.null(characteristics))
    {
        if (!is.mixture.of.normals)
            result <- latentClassMaxDiff(dat, dat$respondent.indices, NULL, n.classes, seed,
                                         initial.parameters, 0, trace, TRUE, lc.tolerance, is.tricked)
        else
            result <- mixtureOfNormalsMaxDiff(dat, n.classes, normal.covariance, seed, initial.parameters,
                                               trace, pool.variance, n.draws, is.tricked)
    }
    else
        result <- varyingCoefficientsMaxDiff(dat, n.classes, seed, initial.parameters, trace, apply.weights,
                                             lc, sub.model.outputs, lc.tolerance, is.tricked)

    result <- accuracyResults(dat, result, questions.left.out)
    if (sub.model.outputs)
        cat("Latent class analysis",
            "BIC:", result$bic,
            "LL:", result$log.likelihood,
            "Number of classes:", n.classes)
    result$n.respondents <- length(dat$respondent.indices) # this should be the unfiltered number of respondents
    result$n.classes <- n.classes
    result$subset <- subset
    result$weights <- weights
    result$n.questions <- dat$n.questions
    result$n.alternatives.per.task <- ncol(dat$X.in)
    result$output <- output
    result$lc <- lc
    result$questions.left.out <- questions.left.out

    resp.pars <- as.matrix(RespondentParameters(result))[dat$subset, ]
    result$respondent.probabilities <- exp(resp.pars) / rowSums(exp(resp.pars))
    result$is.mixture.of.normals <- is.mixture.of.normals
    result$algorithm <- algorithm
    result
}

predictionAccuracies <- function(object, X, n.questions, subset)
{
    score <- rep(NA, nrow(X))
    resp.pars <- as.matrix(RespondentParameters(object)[subset, ])
    n.respondents <- nrow(resp.pars)
    result <- rep(NA, n.respondents)
    for (i in 1:n.respondents)
    {
        pars <- resp.pars[i, ]
        score <- rep(NA, n.questions)
        for (j in 1:n.questions)
        {
            ind <- (i - 1) * n.questions + j
            u <- pars[X[ind, ]]
            score[j] <- all(u[1] > u[-1])
        }
        result[i] <- mean(score)
    }
    result
}

accuracyResults <- function(dat, result, questions.left.out)
{
    n.respondents <- length(dat$respondent.indices)
    in.sample.accuracies <- predictionAccuracies(result, dat$X.in, dat$n.questions.in, dat$subset)
    w <- dat$weights[(1:n.respondents) * dat$n.questions.in]
    result$in.sample.accuracy <- sum(in.sample.accuracies * w) / sum(w)
    if (questions.left.out > 0)
    {
        result$prediction.accuracies <- predictionAccuracies(result, dat$X.out, questions.left.out, dat$subset)
        result$out.sample.accuracy <- sum(result$prediction.accuracies * w) / sum(w)
    }
    else
    {
        result$prediction.accuracies <- in.sample.accuracies
        result$out.sample.accuracy <- NA
    }
    result
}

#' \code{RespondentParameters}
#' @description The parameters for each respondent.
#' @param object A \code{FitMaxDiff} object.
#' @export
RespondentParameters <- function(object)
{
    if (is.null(object$respondent.parameters)) # retained for backwards compatibility
    {
        alternative.names <- if (object$n.classes > 1) rownames(object$coef) else names(object$coef)
        computeRespondentParameters(object, alternative.names, object$input.respondent.pars)
    }
    else
        as.data.frame(object$respondent.parameters)
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
        "MaxDiff: Varying Coefficients"
    else if (x$is.mixture.of.normals)
        "MaxDiff: Mixture of Normals"
    else if (x$algorithm == "HB")
        "MaxDiff: Hierarchical Bayes"
    else
        "MaxDiff: Latent Class Analysis"
    footer <- paste0("n = ", x$n.respondents, "; ")
    if (!is.null(x$subset) && !all(x$subset))
        footer <- paste0(footer, "Filters have been applied; ")
    if (!is.null(x$weights))
        footer <- paste0(footer, "Weights have been applied; Effective sample size: ",
                         FormatWithDecimals(x$effective.sample.size, 2), "; ")
    footer <- paste0(footer, "Number of questions: ", x$n.questions, "; ")
    if (x$questions.left.out > 0)
    {
        footer <- paste0(footer, "Questions used in estimation: ", x$n.questions - x$questions.left.out, "; ")
        footer <- paste0(footer, "Questions left out: ", x$questions.left.out, "; ")
    }
    footer <- paste0(footer, "Alternatives per question: ", x$n.alternatives.per.task, "; ")
    if (x$algorithm != "HB")
    {
        footer <- paste0(footer, "Log-likelihood: ", FormatWithDecimals(x$log.likelihood, 2), "; ")
        footer <- paste0(footer, "BIC: ", FormatWithDecimals(x$bic, 2), "; ")
    }

    footer <- if (!x$lc && x$algorithm != "HB" && !x$is.mixture.of.normals)
        paste0(footer, "Latent class analysis over respondents not applied; ")
    else if (x$is.mixture.of.normals)
    {
        if (x$n.classes == 1)
            paste0(footer, "Mixture of normals: ", x$n.classes, " class; ")
        else
            paste0(footer, "Mixture of normals: ", x$n.classes, " classes; ")
    }
    else if (x$algorithm == "HB")
        footer
    else
    {
        if (x$n.classes == 1)
            paste0(footer, "Latent class analysis: ", x$n.classes, " class; ")
        else
            paste0(footer, "Latent class analysis: ", x$n.classes, " classes; ")
    }

    subtitle <- if (!is.na(x$out.sample.accuracy))
        paste0("Prediction accuracy (leave-", x$questions.left.out , "-out cross-validation): ",
               FormatAsPercent(x$out.sample.accuracy, 3))
    else
        paste0("Prediction accuracy (in-sample): ", FormatAsPercent(x$in.sample.accuracy, 3))

    if (x$n.classes == 1 && is.null(x$covariates.notes)
        && ((!x$is.mixture.of.normals && x$algorithm != "HB") || x$output == "Classes"))
    {
        col.labels <- "Probabilities (%)"
        MaxDiffTableClasses(as.matrix(x$class.preference.shares), col.labels, title, subtitle, footer)
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
        MaxDiffTableClasses(as.matrix(x$class.preference.shares), col.labels, title, subtitle, footer)
    }
}
