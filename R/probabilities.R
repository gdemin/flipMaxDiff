#' \code{dBestWorst}
#' @description Computes the probabilities that a pair of items are chosen as 'best' and 'worst' from a set.
#' @param e.u A vector of \code{exp(utilities)}, where items are shown in the order chosen, from 'best' to 'worst'.
#' @export
dBestWorst <- function(e.u)
{
    k <- length(e.u)
    if (k == 1)
        return(1)
    if (k == 2)
        return(e.u[1] / sum(e.u))
    permutations <- Permutations(2:(k-1))
    possible.rankings <- cbind(1, permutations, k) # The comuputaiton can be made more efficient by computing the first choice separately.
    probs <- apply(possible.rankings, 1, function(x) dExplodedLogit(e.u[x]))
    sum(probs)
}



#' \code{dExplodedLogit}
#' @description Computes the probability of an ordering, where the items are ordered from 'best' to 'worst'.
#' @param e.u A vector of exp(utilities), ordered from most 'preferred to leas'best' to 'worst'.
#' @export
dExplodedLogit <- function(e.u)
{
    prob.best <- e.u[1] / sum(e.u)
    if (length(e.u) == 2)
        return(prob.best)
    prob.best * dExplodedLogit(e.u[-1])
}

#' \code{permutations}
#' @description Computes all possible permitatuions of \code{x}.
#' @param x A vector.
#' @param prefix The empty case; used for recursion - don't play with this.
#' @references http://stackoverflow.com/a/29023189/1547926
Permutations <- function(x, prefix = c() )
{
    if(length(x) == 0)
        return(prefix)
    do.call(rbind, sapply(1:length(x), FUN = function(idx) Permutations( x[-idx], c(prefix, x[idx])), simplify = FALSE))
 }

