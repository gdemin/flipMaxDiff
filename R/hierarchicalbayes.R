#' @importFrom rstan rstan_options stan extract sampling
hierarchicalBayesMaxDiff <- function(dat, n.iterations = 100, n.chains = 1, max.tree.depth = 10,
                                     adapt.delta = 0.8, is.tricked = FALSE)
{
    # We want to replace this call with a proper integration of rstan into this package
    require(rstan)

    # allows Stan chains to run in parallel on multiprocessor machines
    options(mc.cores = parallel::detectCores())

    n.choices <- ncol(dat$X.in)
    n.alternatives <- dat$n.alternatives
    n.respondents <- dat$n
    n.questions <- dat$n.questions
    n.questions.left.in <- dat$n.questions.in

    X <- array(dim = c(n.respondents, n.questions.left.in, n.choices, n.alternatives))
    Y.best <- array(1, dim = c(n.respondents, n.questions.left.in))
    Y.worst <- array(n.choices, dim = c(n.respondents, n.questions.left.in))
    Z <- array(1, dim = c(1, n.respondents))

    for (n in 1:n.respondents)
    {
        for (q in 1:n.questions.left.in)
        {
            X[n, q, , ] <- matrix(0, nrow = n.choices, ncol = n.alternatives)
            for (i in 1:n.choices)
                X[n, q, i, dat$X.in[(n - 1) * n.questions.left.in + q, i]] <- 1
        }
    }

    stan.dat <- list(C = n.choices,
                     K = n.alternatives,
                     R = n.respondents,
                     S = n.questions.left.in,
                     G = 1,
                     YB = Y.best,
                     YW = Y.worst,
                     X = X,
                     Z = Z)

    if (.Platform$OS.type == "unix")
    {
        # Loads a precompiled stan model called mod from sysdata.rda to avoid recompiling.
        # The R code used to generate mod on a linux machine is:
        # mod <- rstan::stan_model(model_code = model.code)
        # devtools::use_data(mod, internal = TRUE, overwrite = TRUE)
        # where model.code is the stan code as a string.
        # Ideally we would want to recompile when the package is built (similar to Rcpp)
        stan.fit <- sampling(mod, data = stan.dat, chains = n.chains, iter = n.iterations)
        # Replace stanmodel with a dummy as stanmodel makes the output many times larger,
        # and is not required for diagnostic plots.
        dummy.stanmodel <- ""
        class(dummy.stanmodel) <- "stanmodel"
        stan.fit@stanmodel <- dummy.stanmodel
    }
    else # windows
    {
        rstan_options(auto_write=TRUE) # writes a compiled Stan program to the disk to avoid recompiling
        stan.fit <- stan(file = "exec/hb.stan", data = stan.dat, iter = n.iterations,
                         chains = n.chains, control = list(max_treedepth = max.tree.depth,
                                                           adapt_delta = adapt.delta))
    }

    resp.pars <- t(colMeans(extract(stan.fit, pars=c("Beta"))$Beta, dims = 1))
    colnames(resp.pars) <- dat$alternative.names

    result <- list(respondent.parameters = resp.pars, stan.fit = stan.fit)
    class(result) <- "FitMaxDiff"
    result
}
