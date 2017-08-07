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

    suppressWarnings(stan.fit <- sampling(mod, data = stan.dat, chains = n.chains, iter = n.iterations))

    resp.pars <- t(colMeans(extract(stan.fit, pars=c("Beta"))$Beta, dims = 1))
    colnames(resp.pars) <- dat$alternative.names

    result <- list(respondent.parameters = resp.pars, stan.fit = stan.fit)
    class(result) <- "FitMaxDiff"
    result
}

.onLoad <- function(libname, pkgname) {

    model.code <- "data {
        int<lower=2> C; // Number of choices in each scenario
        int<lower=1> K; // Number of alternatives
        int<lower=1> R; // Number of respondents
        int<lower=1> S; // Number of scenarios per respondent
        int<lower=0> G; // Number of respondent covariates
        int<lower=1,upper=C> YB[R, S]; // best choices
        int<lower=1,upper=C> YW[R, S]; // worst choices
        matrix[C, K] X[R, S]; // matrix of attributes for each obs
        matrix[G, R] Z; // vector of covariates for each respondent
    }

    parameters {
        matrix[K, R] Beta;
        matrix[K, G] Theta;
        corr_matrix[K] Omega;
        vector<lower=0>[K] tau;
    }
    transformed parameters {
        cov_matrix[K] Sigma = quad_form_diag(Omega, tau);
    }
    model {
        //priors
        to_vector(Theta) ~ normal(0, 10);
        tau ~ cauchy(0, 2.5);
        Omega ~ lkj_corr(2);
        //likelihood
        for (r in 1:R) {
            Beta[,r] ~ multi_normal(Theta*Z[,r], Sigma);
            for (s in 1:S) {
                YB[r,s] ~ categorical_logit(X[r,s] * Beta[,r]);
                YW[r,s] ~ categorical_logit(-X[r,s] * Beta[,r]);
            }
        }
    }"

    message("starting")
    mod <- rstan::stan_model(model_code = model.code)
    devtools::use_data(mod, internal = TRUE, overwrite = TRUE)
    message("finished")
}
