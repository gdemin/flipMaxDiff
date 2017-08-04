varyingCoefficientsMaxDiff <- function(dat, n.classes, seed, initial.parameters, trace, apply.weights, lc,
                                       sub.model.outputs, lc.tolerance, is.tricked)
{
    n.respondents <- length(dat$respondent.indices)
    n.questions.in <- dat$n.questions.in
    alternative.names <- dat$alternative.names
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
            best.n.classes <- NULL
            # Loop over all possible class sizes
            for (n.c in 1:n.levels)
            {
                solution <- latentClassMaxDiff(dat, ind.levels, resp.pars, n.c, seed,
                                               initial.parameters, n.previous.parameters, trace,
                                               apply.weights = apply.weights, lc.tolerance, is.tricked)
                if (solution$bic < best.bic)
                {
                    best.bic <- solution$bic
                    best.solution <- solution
                    best.n.classes <- n.c
                }
            }
            if (best.n.classes > 1)
            {
                resp.pars <- as.matrix(RespondentParameters(best.solution))[dat$subset, ]
                n.previous.parameters <- best.solution$n.parameters
                covariates.notes[i] <- paste0(names(characteristics)[i], " - ", best.n.classes, " classes")
                covariates.chosen <- TRUE
            }
            else
            {
                covariates.notes[i] <- paste0(names(characteristics)[i], " - no selection")
                if (sub.model.outputs)
                    cat("Covariate:", names(characteristics)[i],
                        "BIC:", best.solution$bic,
                        "LL:", best.solution$log.likelihood,
                        "Number of classes:", best.n.classes,
                        (if (best.solution$n.classes == 1) "Excluded" else "Included in final model"),
                        "\n")
            }
        }
    }

    result <- if (lc || is.null(characteristics))
        latentClassMaxDiff(dat, dat$respondent.indices, resp.pars, n.classes, seed, initial.parameters,
                           n.previous.parameters, trace, apply.weights, lc.tolerance, is.tricked)
    else if (covariates.chosen)
        best.solution
    else
        stop("No model applied. Choose different covariates or enable latent class analysis over respondents.")

    result$covariates.notes <- covariates.notes
    result
}
