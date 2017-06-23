context("mixtures of normals")

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("Full covariance", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, n.classes = 2,
                         alternative.names = names, is.mixture.of.normals = TRUE, normal.covariance = "Full")
    expect_error(print(result), NA)
})

test_that("Full covariance with subset and weight", {
    sub <- unclass(tech.data$Q2) <= 3
    wgt <- tech.data$RESPNUM
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, n.classes = 2,
                         alternative.names = names, is.mixture.of.normals = TRUE, normal.covariance = "Full",
                         subset = sub, weights = wgt)
    expect_error(print(result), NA)
})

test_that("Spherical", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, n.classes = 2,
                         alternative.names = names, is.mixture.of.normals = TRUE, normal.covariance = "Spherical")
    expect_error(print(result), NA)
})

test_that("Diagonal", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, n.classes = 2,
                         alternative.names = names, is.mixture.of.normals = TRUE, normal.covariance = "Diagonal")
    expect_error(print(result), NA)
})

test_that("Cross-validation", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, n.classes = 2,
                         alternative.names = names, is.mixture.of.normals = TRUE, normal.covariance = "Full",
                         tasks.left.out = 2)
    expect_error(print(result), NA)
})
