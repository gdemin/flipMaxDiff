context("max-diff logit")

tech.data <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
tech.design <- read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
best <- tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
worst <- tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")

des <- matrix(names[as.integer(as.matrix(tech.design[, -1:-2]))], ncol = 5)

test_that("Estimating logit parameters", {
    # Aggregate
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names)
    q.solution <- c(0, -0.09059317085259,-1.022392901667,0.3712129560247,-0.6596444467744,0.02930605087281,-0.08321462352542,-0.870332546743,-1.105990654593,-0.8674738151854)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    q.probs <- c(0.1358834315826, 0.1241144627641, 0.0488817685569, 0.1969619087625,
                 0.07025650843667, 0.1399245639353, 0.125033634106, 0.05690964591415,
                 0.04496150786246, 0.05707256807926)
    expect_equal(as.vector(result$class.preference.shares), q.probs, tolerance = 0.0001)
    expect_error(print(result), NA)
    # Subset
    sub <- unclass(tech.data$Q2) <= 3
    result = FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, subset = sub, alternative.names = names)
    q.solution <- c(0,-0.09164874531048,-1.416097197804,0.5413188029853,-0.8428037706912,-0.5743677293267,-0.4820668863369,-0.8572557181654,-1.530467590026,-1.143659145101)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.0001)
    # Subset and weight
    wgt <- tech.data$RESPNUM
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, subset = sub, weights = wgt, alternative.names = names)
    q.solution <- c(0,-0.1358158415588,-1.310226780669,0.4291539729508,-0.7617988208062,-0.5332417439411,-0.4662730731195,-0.8116530486171,-1.356435800804,-1.148628546175)
    expect_equal(as.vector(result$coef), q.solution, tolerance = 0.00001)
    expect_equal(result$log.likelihood, -804.218, tolerance = 0.00001)

    expect_error(print(result), NA)

    # Cross validation
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, tasks.left.out = 2)
    expect_error(print(result), NA)
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, subset = sub, weights = wgt, alternative.names = names, tasks.left.out = 2)
    expect_error(print(result), NA)
})

test_that("Latent class", {
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 2)
    q.solution <- structure(c(0, -2.30575675367483, -3.4390021812939, -1.718556454195,
                              -3.03608376429065, -2.43609639597502, -2.42392249235851, -3.28327816635577,
                              -3.50897927456835, -3.44485235322318, 0, 1.38580057447584, 0.428462438590491,
                              1.84870129481078, 0.931693897342025, 1.97524268771123, 1.63144214440238,
                              0.701477372305963, 0.327266346090314, 0.941038167774774), .Dim = c(10L, 2L))
    expect_equal(unname(result$coef), q.solution, tolerance = 0.01)
    expect_equal(result$log.likelihood, -4482.242, tolerance = 0.00001)
    q.probs <- structure(c(0.6118, 0.061, 0.0196, 0.1097, 0.0294, 0.0535, 0.0542,
                           0.0229, 0.0183, 0.0195, 0.0297, 0.1186, 0.0455, 0.1884, 0.0753,
                           0.2138, 0.1516, 0.0598, 0.0412, 0.076), .Dim = c(10L, 2L))
    expect_equal(unname(result$class.preference.shares[, 1:2]), q.probs, tolerance = 0.001)
    expect_error(print(result), NA)

    # Subset and weight
    sub <- unclass(tech.data$Q2) <= 3
    wgt <- tech.data$RESPNUM
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 2, weight = wgt, subset = sub)

    q.solution <- structure(c(0, -1.34109808073123, -3.72839178922714, -0.524605012308453,
                -2.73912325516393, -1.41606753418318, -2.45802001427938, -3.18998838616667,
                -3.78237027293449, -3.3674225091177, 0, 0.517444911638947, 0.275768946657516,
                1.04463763172229, 0.425653550385925, -0.0960750355595251, 0.904803630723954,
                0.858383432709159, 0.19099885898787, 0.164046805635207), .Dim = c(10L, 2L))
    expect_equal(unname(result$coef), q.solution, tolerance = 0.01) # this is quite lax right now. Fit will improve with the addition of a numerical optimization stage
    expect_equal(result$log.likelihood, -730.724, tolerance = 0.01)

    expect_error(print(result), NA)

    expect_error(print(FitMaxDiff(design = tech.design,
                                  version = rep(1, nrow(best)),
                                  best = best,
                                  worst = worst,
                                  alternative.names = names,
                                  n.classes = 2,
                                  output = "Classes")), NA)

    # Cross validation
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 2, tasks.left.out = 3)
    expect_error(print(result), NA)
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 2, weight = wgt, subset = sub, tasks.left.out = 3)
    expect_error(print(result), NA)
})

test_that("Checking some of the inputs", {
    # No Version column in design
    expect_error(FitMaxDiff(design = tech.design[, -1], version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names), NA)
    # No Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -2], version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names), NA)
    # Neither version nor Task column in design
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names), NA)
    # No version
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], best = best, worst = worst, alternative.names = names), NA)
    # inconsistent version information
    expect_error(FitMaxDiff(design = tech.design[, -1:-2], version = 1:nrow(best), best = best, worst = worst, alternative.names = names))
    des <- tech.design
    des$Version <- 3
    expect_error(FitMaxDiff(design = des, best = best, worst = worst, alternative.names = names))
    # No names
    expect_error(FitMaxDiff(design = tech.design, best = best, worst = worst))
    # Incorrect names as a string
    expect_error(suppressWarnings(FitMaxDiff(design = tech.design, best = best, worst = worst, alternative.names = "A,B,C,D,E,F,G,h,I,J")))
    # Correct names as a string
    nms = paste(names, collapse = ", ")
    expect_error(FitMaxDiff(design = tech.design, best = best, worst = worst, alternative.names = nms), NA)
})

test_that("Varying coefficients", {
    # Aggregate
    ll.aggregate <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = NULL, n.classes = 1)$log.likelihood
    # Gender varying coefficient
    expect_error(FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q1), n.classes = 1, lc = FALSE))
    # Apple varying coefficient
    ll.apple <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q3_01), n.classes = 1, lc = FALSE)$log.likelihood
    expect_true(ll.apple > ll.aggregate)
    # Apple varying coefficient with boosting
    ll.apple.boosting.5.classes <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q3_01), n.classes = 5)$log.likelihood
    ll.5.classes <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 5)$log.likelihood
    expect_true(ll.apple.boosting.5.classes > ll.5.classes)

    ll.apple.boosting.1.class <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q3_01), n.classes = 1)$log.likelihood
    expect_true(ll.apple.boosting.1.class > ll.apple)
    expect_true(ll.apple.boosting.1.class > ll.aggregate)

    # Cross validation
    apple.boosting.2.class <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q3_01), n.classes = 2, tasks.left.out = 1)
    expect_error(print(apple.boosting.2.class), NA)

    #Subset
    sub <- unclass(tech.data$Q2) <= 3
    result <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, characteristics = data.frame(tech.data$Q3_01), n.classes = 2, tasks.left.out = 4, subset = sub)
    expect_error(print(result), NA)
})


test_that("Saving variables", {
    sub <- c(FALSE, FALSE, rep(TRUE, 100), rep(FALSE, 200))
    # Posterior probabilities.
    lc.3 <- FitMaxDiff(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst, alternative.names = names, n.classes = 3, subset = sub)
    pp <- lc.3$posterior.probabilities[sub, ]
    expect_equal(ncol(pp), 3)
    expect_equal(sd(apply(pp, 1, sum)), 0)
    # Individual-level paramters
    pars <- RespondentParameters(lc.3)
    expect_equal(sum(!is.na(pars[, 1])), sum(sub))
    # Segment memberships
    m <- table(Memberships(lc.3))
    expect_equal(length(m), 3)
    expect_equal(sum(m), sum(sub))
})

test_that("Experimental designs with versions", {
    # President.
    dat <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/61/President.sav", to.data.frame = TRUE))
    des <- read.csv("http://wiki.q-researchsoftware.com/images/9/9d/PresidentialDesign.csv")
    best <- dat[, c("MDmost_1", "MDmost_2", "MDmost_3", "MDmost_4"   ,"MDmost_5"  , "MDmost_6", "MDmost_7"  ,"MDmost_8","MDmost_9","MDmost_10" )]
    worst <- dat[, c("MDleast_1", "MDleast_2", "MDleast_3" , "MDleast_4",  "MDleast_5"  ,"MDleast_6",  "MDleast_7", "MDleast_8",  "MDleast_9", "MDleast_10")]
    names <- c("Decent/ethical", "Plain-speaking", "Healthy", "Successful in business", "Good in a crisis", "Experienced in government", "Concerned for minorities", "Understands economics", "Concerned about global warming", "Concerned about poverty", "Has served in the military", "Multilingual", "Entertaining", "Male", "From a traditional American background", "Christian")
    expect_error(suppressWarnings(FitMaxDiff(design = des, dat$MDversion, best = best, worst = worst, alternative.names = names)))
    names <- c("Decent/ethical", "Plain-speaking", "Healthy", "Successful in business", "Good in a crisis", "Experienced in government", "Concerned for the welfare of minorities", "Understands economics", "Concerned about global warming", "Concerned about poverty", "Has served in the military", "Multilingual", "Entertaining", "Male", "From a traditional American background", "Christian")
    expect_error(FitMaxDiff(design = des, dat$MDversion, best = best, worst = worst, alternative.names = names), NA)
    expect_error(FitMaxDiff(design = des, dat$MDversion, best = best, worst = worst, alternative.names = names, n.classes = 3), NA)
    # Example from Q wiki
    # dat <- suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/6/66/MaxDiffSetupExample.sav", to.data.frame = TRUE))
    # des <- read.csv("http://wiki.q-researchsoftware.com/images/2/24/ExampleMaxDiffDesign.csv")
    # best <- dat[, paste0("Q11task", 1:13,"most")]
    # worst <- dat[, paste0("Q11task", 1:13,"least")]
    # names <- paste("Brand", 1:13)
    # expect_error(FitMaxDiff(design = des, dat$MDversion, best = best, worst = worst, alternative.names = names))

})


