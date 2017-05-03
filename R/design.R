#' \code{MaxDiffDesign}
#' @description Creates an experimental design for a max-diff experiment.
#' @param number.alternatives The number of alternatives in the experiment. For example, if you are doing a study investigating preferences for 10 brands, then 10 is the number of alternatives.
#' @param number.questions The number of max-diff questions to show to respondents. Sawtooth Software suggests that a rough guideline is: \code{Number of questions >= 3 * Number of alternatives / Alternatives per question}.
#' @param alternatives.per.question For example, if you have a study of 10 brands, and in each question you show five brands, asking the respondent to choose the one of the five that they like the most and the one that they like the least, then \code{Alternatives per question = 5}.
#' @param n.repeats The number of times that the algorithm seeks to find a solution. The higher the number, the greater the chance that the best possible solution is found. For most problems, this makes little difference (i.e., a marginally sub-optimal experimental design will tend not to have any meaningful consequence on the conclusions drawn from the analyses).
#' @import AlgDesign
#' @param x The experimental design for a respondent  (a \code{\link{list}}).
#' @export
MaxDiffDesign <- function(number.alternatives, number.questions, alternatives.per.question, n.repeats = 1000){
    # Check that the parameters are appropriate
    # Sawtooth recommends that number.questions >= 3 * number.alternatives / alternatives.per.question
    if (number.questions < 3 * number.alternatives / alternatives.per.question)
        warning("It is recomended that number.questions >= 3 * number.alternatives / alternatives.per.question");
    best.result = NULL
    best.D = -Inf
    for (i in 1:n.repeats){
        alg.results <- optBlock(~.,withinData=factor(1:number.alternatives),
                                blocksizes=rep(alternatives.per.question,number.questions),
                                nRepeats=5000) #BIB
        if (alg.results$D > best.D){
            best.result = alg.results
            best.D = alg.results$D
        }
    }
    design <- matrix(best.result$rows, nrow = number.questions, byrow = TRUE, dimnames = list(Question = 1:number.questions, Alternative = 1:alternatives.per.question))
    CheckMaxDiffDesign(design)
}


#' \code{CheckMaxDiffDesign}
#' @description Produces summary statistics for a max-diff design.
#' @param design A \link{\code{matrix}}, where each row represents a question or task, and each column
#' shows the alternatives to be shown.
#' @export
CheckMaxDiffDesign <- function(design)
{
    design <- as.matrix(design)
    number.questions <- nrow(design)
    number.alternatives <- max(design)
    alternatives.per.question <- ncol(design)
    binary.design <- matrix(0,number.questions,number.alternatives, dimnames = list(Question = paste("Question", 1:number.questions), Alternative = 1:number.alternatives))
    for (q in 1:number.questions)
        binary.design[q, design[q, ]] <- 1
    n.appearances.per.alternative <- table(as.numeric(design))
    combinations.of.alternatives <- crossprod(table(c(rep(1:number.questions, rep(alternatives.per.question,number.questions))), as.integer(t(design))))
    colnames(binary.design) <- c("Alternative 1", 2:number.alternatives)
    list(binary.design = binary.design,
         design = design,
         frequencies = n.appearances.per.alternative,
         pairwise.frequencies=combinations.of.alternatives,
         binary.correlations = round(cor(binary.design),2))
}
