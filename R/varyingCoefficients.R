varyingCoefficientsMaxDiff <- function(dat, alternative.names, n.classes, seed, initial.parameters,
                                       trace, apply.weights, lc, sub.model.outputs)
{
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

    result$covariates.notes <- covariates.notes
    result
}
