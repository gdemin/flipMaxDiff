#' @importFrom rstan rstan_options stan extract
hierarchicalBayesMaxDiff <- function(dat, n.iterations = 100, n.chains = 1, is.tricked = FALSE)
{
    # We want to replace this call with a proper integration of rstan into this package
    require(rstan)

    # writes a compiled Stan program to the disk to avoid recompiling
    rstan_options(auto_write=TRUE)
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

    suppressWarnings(stan.fit <- stan(file = "exec/hb.stan",
                                      data = stan.dat,
                                      iter = n.iterations,
                                      chains = n.chains))

    resp.pars <- t(colMeans(extract(stan.fit, pars=c("Beta"))$Beta, dims = 1))
    colnames(resp.pars) <- dat$alternative.names

    result <- list(respondent.parameters = resp.pars, stan.fit = stan.fit)
    class(result) <- "FitMaxDiff"
    result
}
