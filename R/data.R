#' \code{IntegrateDesignAndData}
#' @description A matrix where in each row, the first element is the 'best' choice and the last element is the
#' 'worst' choice.
#' @param design A \code{data.frame}, where the first variable is called 'Version', the second is called 'Task',
#' and the remaining variables contain the alternatives shown in each tas.
#' @param version A vector of integers showing the version of the design shown to each respondent.
#' @param best A matrix of integers showing the choices made by each respondent on each of the tasks. One column
#' for each task. The integers need to corresponde to the \code{design} vector of integers showing the version of
#' the design shown to each respondent. Coerced to a matrix if a \code{data.frame}.
#' @param worst  A matrix of integers showing the choice of 'worst'.
#' @export
IntegrateDesignAndData <- function(design, version, best, worst)
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
            respondent.indices[[i]] <- (1:n.tasks) + (i - 1) * n.tasks
            c = c + 1
        }
    }
    list(X = X, respondent.indices = respondent.indices)
}
