context("design")
test_that("design",
{
    des <- MaxDiffDesign(number.alternatives = 10,
                            number.questions = 6,
                            alternatives.per.question = 5,
                            n.repeats = 1)
    tech.design = read.csv("http://wiki.q-researchsoftware.com/images/7/78/Technology_MaxDiff_Design.csv")
    z1 <- CheckMaxDiffDesign(tech.design[, -1:-2])$binary.correlations
    z1 <- range(z1[lower.tri(z1)])
    z2 <- des$binary.correlations
    z2 <- range(z2[lower.tri(z2)])
    expect_equal(z1, z2)
    # warnings
    expect_warning(MaxDiffDesign(number.alternatives = 10,
                            number.questions = 6,
                            alternatives.per.question = 5,
                            n.repeats = 1), NA)
    expect_warning(MaxDiffDesign(number.alternatives = 10,
                            number.questions = 6,
                            alternatives.per.question = 4,
                            n.repeats = 1))

})

