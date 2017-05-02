context("probabilities")

test_that("Permutations", {
    expect_equal(nrow(Permutations(1:3)), 6)
    expect_equal(nrow(Permutations(1:4)), 24)
    expect_equal(nrow(Permutations(1:5)), 120)
})

test_that("dExplodedLogit", {
    expect_equal(dExplodedLogit(c(1, 1)), 0.5)
    expect_equal(dExplodedLogit(c(1, 1, 1)), 1 / 6)
})

test_that("Single task", {
    ## Example 1
    # Utilities
    zb = 4:1
    # Computing the probability direclty
    p.best = exp(zb[1]) / sum(exp(zb[1:4]))
    p.not.worst = exp(zb[2]) / sum(exp(zb[2:4])) * exp(zb[3]) / sum(exp(zb[3:4])) + exp(zb[3]) / sum(exp(zb[2:4])) * exp(zb[2]) / sum(exp(zb[c(2, 4)]))
    p.direct = p.best * p.not.worst
    # Computing the probabilities after some trivial substitution
    p.sub = exp(4) / sum(exp(1:4)) * (exp(3) / sum(exp(3:1)) * exp(2) / sum(exp(2:1)) + exp(2) / sum(exp(3:1)) * exp(3) / sum(exp(c(3, 1)) ))
    # Computing the probailities using the biased urn
    .probabilityUsingBiasedUrn = function(b){
      eb = exp(b)
      k = length(eb)
      d.best = eb[1]/sum(eb)
      d.not.worst = BiasedUrn::dMWNCHypergeo(c(rep(1,k-2),0), rep(1,k-1),k-2,eb[-1], precision = 1E-7)
      d.best * d.not.worst}
    p.biased.urn = .probabilityUsingBiasedUrn(zb)
    # Computing probabilities exactly using simple function
    p.simple = dBestWorst(exp(zb)) # exp as this is how it expects is parameters
    # The tests
    expect_equal(p.direct, p.sub, p.biased.urn, p.simple)

    ## Example 2
    # Utilities
    set.seed(1232)
    zb = runif(4)
    # Computing the probability direclty
    p.best = exp(zb[1]) / sum(exp(zb[1:4]))
    p.not.worst = exp(zb[2]) / sum(exp(zb[2:4])) * exp(zb[3]) / sum(exp(zb[3:4])) + exp(zb[3]) / sum(exp(zb[2:4])) * exp(zb[2]) / sum(exp(zb[c(2, 4)]))
    p.direct = p.best * p.not.worst
    # Computing the probailities using the biased urn
    .probabilityUsingBiasedUrn = function(b){
      eb = exp(b)
      k = length(eb)
      d.best = eb[1]/sum(eb)
      d.not.worst = BiasedUrn::dMWNCHypergeo(c(rep(1,k-2),0), rep(1,k-1),k-2,eb[-1], precision = 1E-7)
      d.best * d.not.worst}
    p.biased.urn = .probabilityUsingBiasedUrn(zb)
    # Computing probabilities exactly using simple function
    p.simple = dBestWorst(exp(zb)) # exp as this is how it expects is parameters
    # The tests
    expect_equal(p.direct, p.biased.urn, p.simple)

    ## Examples > 2
    # Utilities
    set.seed(1232)
    for (i in 2:10)
    {
        zb = runif(i)
        # Computing the probailities using the biased urn
        .probabilityUsingBiasedUrn = function(b){
          eb = exp(b)
          k = length(eb)
          d.best = eb[1]/sum(eb)
          if (k == 2)
              return(d.best)
          d.not.worst = BiasedUrn::dMWNCHypergeo(c(rep(1,k-2),0), rep(1,k-1),k-2,eb[-1], precision = 1E-7)
          d.best * d.not.worst}
        p.biased.urn = .probabilityUsingBiasedUrn(zb)
        # Computing probabilities exactly using simple function
        p.simple = dBestWorst(exp(zb)) # exp as this is how it expects is parameters
        # The tests
        expect_equal(p.biased.urn, p.simple)
    }
})






