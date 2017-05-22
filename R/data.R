#' \code{IntegrateDesignAndData}
#' @description A matrix where in each row, the first element is the 'best' choice and the last element is the
#' 'worst' choice.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst A matrix of integers showing the choice of 'worst'.
#' @param seed Seed for cross validation
#' @param tasks.left.out Number of tasks to leave out for cross-validation.
#' @export
IntegrateDesignAndData <- function(design, version, best, worst, seed, tasks.left.out = 0)
{
    if (is.data.frame(best))
        best = as.matrix(suppressWarnings(flipTransformations::AsNumeric(best, binary = FALSE)))
    if (is.data.frame(worst))
        worst = as.matrix(suppressWarnings(flipTransformations::AsNumeric(worst, binary = FALSE)))
    n <- length(version)
    if (n != nrow(best) | n != nrow(worst))
        stop("'version', 'best', and 'worst', all need to have the same sample size.")
    n.tasks <- ncol(best)
    if (n.tasks != ncol(worst))
        stop("'best' and 'worst' need to have the same number of columns.")
    n.alternatives <- ncol(design) - 2
    if (max(best) != max(worst) | max(best) != n.alternatives)
        stop(paste0("There are ", n.alternatives, " alternatives per task in the design; this should be the maximum value in 'best' and 'worst' (but is not)."))
    X <- matrix(NA, ncol = n.alternatives, nrow = n * n.tasks)
    respondent.indices = vector("list", n)
    c = 1
    for (i in 1:n)
    {
        respondent.design <- as.matrix(design[design$Version == version[i], -1:-2])
        if (nrow(respondent.design) !=  ncol(best))
            stop("The 'design' has a different number of tasks to the 'best' and 'worst' data")
        respondent <- vector("list", n.tasks)
        for (t in 1:n.tasks)
        {
            task.design <- respondent.design[t, ]
            b.position <- best[i, t]
            w.position <- worst[i, t]
            positions <- c(b.position, (1:n.alternatives)[c(-b.position, -w.position)], w.position)
            X[c, ] <- task.design[positions]
            c = c + 1
        }
    }

    if (tasks.left.out == 0)
    {
        X.in <- X
        X.out <- NULL
        for (i in 1:n)
            respondent.indices[[i]] <- (1:n.tasks) + (i - 1) * n.tasks
    }
    else
    {
        left.out <- leftOutTasks(n, n.tasks, tasks.left.out, seed)
        X.in <- X[!left.out, ]
        X.out <- X[left.out, ]
        tasks.left.in <- n.tasks - tasks.left.out
        for (i in 1:n)
            respondent.indices[[i]] <- (1:tasks.left.in) + (i - 1) * tasks.left.in
    }

    list(X.in = X.in, X.out = X.out, respondent.indices = respondent.indices)
}

#' @importFrom flipData CalibrateWeight CleanSubset CleanWeights
#' @importFrom flipFormat TrimWhitespace
cleanAndCheckData <- function(design, version, best, worst, alternative.names, subset = NULL, weights = NULL,
                              characteristics = NULL, seed, tasks.left.out = 0)
{
    # Check for alternative formats of design, and coerce if not the standard

    # Binary design
    if (class(design) == "matrix")
    {
        if (class(design[1, 1]) == "logical")
        {
            design <- t(apply(design, 1, which))
        } else if (length(which(!design %in% c(0,1))) == 0) {
            design <- design == 1
            design <- t(apply(design, 1, which))
        }
        design <- as.data.frame(design)
        names(design) <- paste0("Alt.", 1:ncol(design))
    } else if (class(design) == "list")
    {
        if (!is.null(design$design))
            design <- as.data.frame(design$design)
        else
            stop("The design should be a data frame.")
    }

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
    if (tasks.left.out >= n.tasks)
        stop("The number of questions left out must be less than the total number of questions.")
    tasks.left.in <- n.tasks - tasks.left.out
    weights <- if (is.null(weights))
        rep(1, tasks.left.in * length(version))
    else
        rep(weights, each = tasks.left.in)
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
         n.tasks = n.tasks,
         n.tasks.in = tasks.left.in,
         respondent.indices = dat$respondent.indices,
         characteristics = characteristics,
         subset = subset)
}

leftOutTasks <- function(n.respondents, n.tasks, n.tasks.left.out, seed)
{
    set.seed(seed)
    sapply(rep(n.tasks, n.respondents), function(x) (1:n.tasks) %in% sample(x, n.tasks.left.out))
}
