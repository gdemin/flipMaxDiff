context("data")

test_that("Reading data works", {

    tech.data = suppressWarnings(foreign::read.spss("http://wiki.q-researchsoftware.com/images/f/f1/Technology_2017.sav", to.data.frame = TRUE))
    tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")

    best = tech.data[, c("Q5a_left", "Q5b_left", "Q5c_left", "Q5d_left", "Q5e_left", "Q5f_left")]
    worst = tech.data[, c("Q5a_right", "Q5b_right", "Q5c_right", "Q5d_right", "Q5e_right", "Q5f_right")]
    # expect_error(IntegrateDesignAndData(design = tech.design, version = rep(1, nrow(best)), best = best, worst = worst), NA)

    names <- c("Apple", "Microsoft", "IBM", "Google", "Intel", "Samsung", "Sony", "Dell", "Yahoo", "Nokia")
    list.design = MaxDiffDesign(number.alternatives = 10, number.questions = 6, alternatives.per.question = 5, n.repeats = 1)
    expect_error(cleanAndCheckData(design = list.design, best = best, worst = worst, alternative.names = names))
    list.design$design <- tech.design
    expect_error(cleanAndCheckData(design = list.design, best = best, worst = worst, alternative.names = names), NA)

    binary.design = list.design$binary.design
    binary.design[binary.design == 1] <- 0
    for (r in 1:nrow(binary.design))
        for (c in 1:5)
            binary.design[r, c] = tech.design[r, c + 2]
    expect_error(cleanAndCheckData(design = binary.design, best = best, worst = worst, alternative.names = names), NA)
})

