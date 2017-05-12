context("max-diff latent class")

tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names = c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("Estimating logit parameters", {
    # weight
    wgt = tech.data$RESPNUM
    result = LatentClassMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 2, weights = wgt)
    q.solution <- c(-2.386, -3.5968, -1.8937, -3.1371, -2.3911, -2.5816, -3.4425, -3.6203, -3.5637, # class 2 parameters in Q
                    1.361, 0.5023, 1.8976, 1.1576, 2.001, 1.7042, 0.8686, 0.3861, 0.8134, # class 1 parameters in Q
                    0.3314001)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)
    expect_equal(result$log.likelihood, -4447.059, tolerance = 0.00001)
})
