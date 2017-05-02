context("max-diff logit")

test_that("Estimating logit parameters", {

    tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
    tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")

    best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
    worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
    names = c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, names = names)
    q.solution <- c(0, -0.09059317085259,-1.022392901667,0.3712129560247,-0.6596444467744,0.02930605087281,-0.08321462352542,-0.870332546743,-1.105990654593,-0.8674738151854)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)

})

