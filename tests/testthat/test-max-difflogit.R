context("max-diff logit")

tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names = c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

test_that("Estimating logit parameters", {
    # Aggregate
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, names = names)
    q.solution <- c(0, -0.09059317085259,-1.022392901667,0.3712129560247,-0.6596444467744,0.02930605087281,-0.08321462352542,-0.870332546743,-1.105990654593,-0.8674738151854)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    # Subset
    sub = unclass(tech.data$Q2) <= 3
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, subset = sub, names = names)
    q.solution <- c(0,-0.09164874531048,-1.416097197804,0.5413188029853,-0.8428037706912,-0.5743677293267,-0.4820668863369,-0.8572557181654,-1.530467590026,-1.143659145101)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    # Subset and weight
    wgt = tech.data$RESPNUM
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, subset = sub, weights = wgt, names = names)
    q.solution <- c(0,-0.1358158415588,-1.310226780669,0.4291539729508,-0.7617988208062,-0.5332417439411,-0.4662730731195,-0.8116530486171,-1.356435800804,-1.148628546175)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)
    expect_equal(result$log.likelihood, -804.218, tolerance = 0.00001)
})

test_that("Checking some of the inputs", {
    # No Version column in design
    expect_error(FitMaxDiff(design = tech.design[, -1], version = rep(1, nrow(best)), best = best, worst = worst, names = names), NA)
    # No Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -2], version = rep(1, nrow(best)), best = best, worst = worst, names = names), NA)
    # Neither version nor Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], version = rep(1, nrow(best)), best = best, worst = worst, names = names), NA)
    # No version
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], best = best, worst = worst, names = names), NA)
    # inconsistent version information
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], version = 1:nrow(best), best = best, worst = worst, names = names))
    des <- tech.design
    des$Version <- 3
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], best = best, worst = worst, names = names))
    # No names
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], best = best, worst = worst), NA)

    q.solution <- c(0,-0.1358158415588,-1.310226780669,0.4291539729508,-0.7617988208062,-0.5332417439411,-0.4662730731195,-0.8116530486171,-1.356435800804,-1.148628546175)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)
    expect_equal(result$log.likelihood, -804.218, tolerance = 0.00001)
})


