#' \code{IntegrateDesignAndData}
#' @description A matrix where in each row, the first element is the 'best' choice and the last element is the
#' 'worst' choice.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Question',
#' and the remaining variables contain the alternatives shown in each question.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best Amatrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to correspond to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst As with 'best', except denoting worst..
#' @param names The names usedof the alterantives.
#' @param seed Seed for cross validation
#' @param tasks.left.out Number of tasks to leave out for cross-validation.
#' @export
IntegrateDesignAndData <- function(design, version, best, worst, seed, tasks.left.out = 0)
{
    n <- length(version)
    if (n != nrow(best) | n != nrow(worst))
        stop("'version', 'best', and 'worst', all need to have the same sample size.")
    n.questions <- ncol(best)
    if (n.questions != ncol(worst))
        stop("'best' and 'worst' need to have the same number of columns.")
    n.alternatives <- ncol(design) - 2
    X <- matrix(NA, ncol = n.alternatives, nrow = n * n.questions)
    respondent.indices = vector("list", n)
    c = 1
    for (i in 1:n)
    {
        respondent.design <- as.matrix(design[design$Version == version[i], -1:-2])
        if (nrow(respondent.design) !=  ncol(best))
            stop("The 'design' has a different number of tasks to the 'best' and 'worst' data")
        respondent <- vector("list", n.questions)
        for (t in 1:n.questions)
        {
            task.design <- respondent.design[t, ]
            b <- best[i, t]
            w <- worst[i, t]
            rnking <- c(b, task.design[!task.design %in% c(b, w)], w)
            if (length(rnking) != length(task.design))
                stop(paste0("There are choices in the 'best' and 'worst' data that are not consistent with the design.
                     Respondent ", i, ", question ", t, "."))
            X[c, ] <- rnking
            c = c + 1
        }
    }

    if (tasks.left.out == 0)
    {
        X.in <- X
        X.out <- NULL
        for (i in 1:n)
            respondent.indices[[i]] <- (1:n.questions) + (i - 1) * n.questions
    }
    else
    {
        left.out <- leftOutTasks(n, n.questions, tasks.left.out, seed)
        X.in <- X[!left.out, ]
        X.out <- X[left.out, ]
        tasks.left.in <- n.questions - tasks.left.out
        for (i in 1:n)
            respondent.indices[[i]] <- (1:tasks.left.in) + (i - 1) * tasks.left.in
    }

    list(X.in = X.in, X.out = X.out, respondent.indices = respondent.indices)
}

#' @importFrom flipData CalibrateWeight CleanSubset CleanWeights
#' @importFrom flipFormat TrimWhitespace
cleanAndCheckData <- function(design, version = NULL, best, worst, alternative.names, subset = NULL, weights = NULL,
                              characteristics = NULL, seed, tasks.left.out = 0)
{
    # Tidying up the best and worst choices.
    if (is.factor(best[[1]]))
    {
        best <- trimws(sapply(best, as.character))
        worst <- trimws(sapply(worst, as.character))
    }
    alternative.names <- trimws(alternative.names)
    if (is.character(best[[1]]))
    {
        best <- apply(best, 2, function(x) match(x, alternative.names))
        worst <- apply(worst, 2, function(x) match(x, alternative.names))
    }
    ## Check for alternative formats of design, and coerce if not the standard
    # Binary design
    if (diff(range(design)) == 1)
    {
        if (class(design[1, 1]) != "logical")
            design <- design == 1
        design <- t(apply(design, 1, which))
        names(design) <- paste0("Alt.", 1:ncol(design))
    }
    else if (class(design) == "list")
        design <- if(is.null(design$versions.design)) design$design else design$versions.design
    design <- as.data.frame(design)

    # Renaming 'Task' as 'Question'
    task.index <- match("Task", names(design))
    if (!is.na(task.index))
        names(design)[task.index] <- "Question"
    # Cleaning and checking data
    n <- nrow(best)
    if (!is.null(design$Version))
    {
        t <- table(design$Version)
        if (any(t != t[1]))
            stop("Versions need to have the same number of tasks.")
        else
            n.questions <- unname(t[1])
    } else
        n.questions <- nrow(design)
    if (!is.null(weights))
        weights <- CleanWeights(weights)
    subset <- CleanSubset(subset, n)
    n <- sum(subset)
    if (!is.null(weights))
    {
        weights <- weights[subset]
        weights <- CalibrateWeight(weights)
    }
    if (tasks.left.out >= n.questions)
        stop("The number of questions left out must be less than the total number of questions.")
    tasks.left.in <- n.questions - tasks.left.out
    weights <- if (is.null(weights))
        rep(1, tasks.left.in * n) else rep(weights, each = tasks.left.in)
    best <- as.data.frame(best[subset, ])
    worst <- as.data.frame(worst[subset, ])
    if (is.null(design$Version))
        design <- cbind(Version = 1, design)
    if (is.null(version))
        version <- rep(1, n)
    version <- version[subset]
    versions.in.variable <- sort(unique(version))
    versions.in.design <- sort(unique(design$Version))
    compare.version <- all.equal(versions.in.variable, versions.in.design)
    if (is.character(compare.version))
        stop("The 'design' and 'version' have incompatible version numbers.")
    if (is.null(design$Question))
        design <- cbind(Question = 1:n.questions, design)
    dat <- IntegrateDesignAndData(design, version, best, worst, seed, tasks.left.out)
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

    list(X.in = dat$X.in,
         X.out = dat$X.out,
         weights = weights,
         alternative.names = alternative.names,
         n = n,
         n.alternatives = n.alternatives,
         n.questions = n.questions,
         n.questions.in = tasks.left.in,
         respondent.indices = dat$respondent.indices,
         characteristics = characteristics,
         subset = subset)
}

leftOutTasks <- function(n.respondents, n.questions, n.questions.left.out, seed)
{
    set.seed(seed)
    sapply(rep(n.questions, n.respondents), function(x) (1:n.questions) %in% sample(x, n.questions.left.out))
}
