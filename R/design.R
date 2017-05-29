#' \code{MaxDiffDesign}
#' @description Creates an experimental design for a max-diff experiment.
#' @param number.alternatives The number of alternatives in the experiment. For example, if you are doing a study investigating preferences for 10 brands, then 10 is the number of alternatives.
#' @param number.questions The number of max-diff questions to show to respondents. Sawtooth Software suggests that a rough guideline is: \code{Number of questions >= 3 * Number of alternatives / Alternatives per question}.
#' @param alternatives.per.question For example, if you have a study of 10 brands, and in each question you show five brands, asking the respondent to choose the one of the five that they like the most and the one that they like the least, then \code{Alternatives per question = 5}. That is, the number of options shown in each question.
#' @param n.repeats The number of times that the algorithm seeks to find a solution. The higher the number, the greater the chance that the best possible solution is found. For most problems, this makes little difference (i.e., a marginally sub-optimal experimental design will tend not to have any meaningful consequence on the conclusions drawn from the analyses).
#' @param n.versions The number of versions of the experimental design (defaults to 1). Subsequent versions are obtained by permutting the columns of the binary design.
#' @param seed Random number seed for generation of the experimental design.
#' @import AlgDesign
#' @export
MaxDiffDesign <- function(number.alternatives, number.questions, alternatives.per.question, n.versions = 1, n.repeats = 1000, seed = 1223){
    # Check that the parameters are appropriate
    # Sawtooth recommends that number.questions >= 3 * number.alternatives / alternatives.per.question
    if (alternatives.per.question >= number.alternatives)
        stop("The number of alternatives per question must be less than the number of alternatives.")
    set.seed(seed)
    best.result <- NULL
    best.D <- -Inf
    for (i in 1:n.repeats){

        alg.results <- try(optBlock(~.,withinData=factor(1:number.alternatives),
                                    blocksizes=rep(alternatives.per.question,number.questions),
                                    nRepeats=5000), silent = TRUE) #BIB, silent = TRUE))
        if (any("try-error" %in% class(alg.results)))
            stop("Unable to compute experimental design. It is likely that your inputs are not sensible.")
        if (alg.results$D > best.D)
        {
            best.result = alg.results
            best.D = alg.results$D
        }
    }
    design <- matrix(best.result$rows, nrow = number.questions, byrow = TRUE, dimnames = list(Questions = paste("Question", 1:number.questions), Alternatives = paste("Option", 1:alternatives.per.question)))
    result <- CheckMaxDiffDesign(design)
    if (n.versions > 1)
        result <- c(result, multipleVersionDesign(result, n.versions))
    result
}


#' \code{multipleVersionDesign}
#' @param original The experiment design that is randomized to individual-level vesions.
#' @param n.versions The number of versions of the experimental design (defaults to 1). Subsequent versions are obtained by permutting the columns of the binary design.
multipleVersionDesign <- function(original, n.versions)
{
    binary.design <- original$binary.design
    number.alternatives <- ncol(binary.design)
    # Creating an array to store the outputs
    if (is.matrix(original))
       stop("Select 'Detailed outputs' on the experimental design.")
    alternatives.per.question <- ncol(original$design)
    number.questions <- nrow(original$design)
    nrows <- number.questions * n.versions
    cnames <- colnames(original$design)
    randomized.designs <- matrix(NA,
         nrow = nrows, ncol = 2 + alternatives.per.question,
         dimnames = list(1:nrows, c("Version", "Question", cnames)))
    randomized.designs[, 1] <- rep(1:n.versions, each = number.questions)
    randomized.designs[, 2] <- rep(1:number.questions)
    randomized.designs[1:number.questions, -1:-2] <- original$design
    big.binary.design <- matrix(NA, nrows, number.alternatives, dimnames = list(1:nrows, colnames(binary.design)))
    big.binary.design[1:number.questions, ] <- binary.design
    # Randomly rearranging the columns.
    set.seed(1223)
    for (i in 2:n.versions)
    {
        rows <- (i - 1) * number.questions + 1:number.questions
        d <- binary.design[, sample(1:number.alternatives, number.alternatives, replace = FALSE)]
        big.binary.design[rows, ] <- d
        d <- t(d)
        design <- matrix(row(d)[d == 1], byrow = TRUE, ncol = alternatives.per.question)
        randomized.designs[rows, -1:-2] <- design
    }
    # Summmary statistics
    correlations <- round(cor(big.binary.design), 2)
    pairwise.frequencies <- crossprod(big.binary.design)
    dimnames(pairwise.frequencies) <- dimnames(correlations) <- list(Alternative = 1:number.alternatives, Alternative = 1:number.alternatives)
    list(versions.binary.correlations = correlations,
         versions.pairwise.frequencies = pairwise.frequencies,
         versions.design = randomized.designs)
}


#' \code{CheckMaxDiffDesign}
#' @description Produces summary statistics for a max-diff design.
#' @param design A \code{\link{matrix}}, where each row represents a question or task, and each column
#' shows the alternatives to be shown.
#' @export
CheckMaxDiffDesign <- function(design)
{
    design <- as.matrix(design)
    number.questions <- nrow(design)
    number.alternatives <- max(design)
    alternatives.per.question <- ncol(design)
    binary.design <- matrix(0,number.questions,number.alternatives, dimnames = list(Question = paste("Question", 1:number.questions), Alternative = 1:number.alternatives))
    colnames(binary.design) <- c("Alternative 1", 2:number.alternatives)
    if (number.questions < 3 * number.alternatives / alternatives.per.question)
        warning(paste0("You have specified ", number.questions, " questions. It is sometimes recommended that number.questions >= 3 * number.alternatives / alternatives.per.question (i.e., that you should have at least ", ceiling(3 * number.alternatives / alternatives.per.question), " questions)."))
    for (q in 1:number.questions)
        binary.design[q, design[q, ]] <- 1
    n.appearances.per.alternative <- table(as.numeric(design))
    if ((min.a <- min(n.appearances.per.alternative)) < 3)
        warning(paste0("One or more of your alternatives appears only ", min.a, " two times. A common recommendation is that each alternative should appear 3 times. You can review the frequencies by viewing the detailed outputs."))
    if (min.a != max(n.appearances.per.alternative))
        warning(paste0("Your design is not balanced. That is, some alternatives appear more frequently than others. You can review the frequencies by viewing the detailed outputs."))
    correlations <- round(cor(binary.design), 2)
    cors <- abs(correlations[lower.tri(correlations)])
    cor.max  <- max(cors, na.rm = TRUE)
    cor.min <- min(cors, na.rm = TRUE)
    if (any(is.na(cors)))
        warning("Some of the binary correlations are zero. This is only a problem of the cause is not that an alternative always appears in the design.")
    if (cor.max > 0.5)
        warning(paste0("The largest binary absolute correlation is ", cor.max, ". You should consider having more questions. You can review the binary correlations by viewing the detailed outputs."))
    if (cor.min != cor.max)
        warning(paste0("The absolute value of the correlations varies from ", cor.min, " to ", cor.max, ". This may not be a problem, but ideally the absolute value of the correlations should be constant (this is not always possible). Consider increasing the number of questions."))
    pairwise.frequencies <- crossprod(binary.design)
    min.pairwise <- min(pairwise.frequencies)
    if (min.pairwise == 0)
        warning(paste0("Some alternatives never appear together. You can review the pairwise frequencies by viewing the detailed outputs."))
    appearance.ratios <- sweep(pairwise.frequencies, 1, n.appearances.per.alternative, "/")
    if (any(appearance.ratios[lower.tri(appearance.ratios)] == 1))
        warning(paste0("Some alternatives only ever appear together. You can review the pairwise frequencies by viewing the detailed outputs."))
    dimnames(pairwise.frequencies) <- dimnames(correlations) <- list(Alternative = 1:number.alternatives, Alternative = 1:number.alternatives)
    list(binary.design = binary.design,
         design = design,
         frequencies = n.appearances.per.alternative,
         pairwise.frequencies = pairwise.frequencies,
         binary.correlations = correlations)
}
