#' Adjust p-values to control FWER given permutations (Westfall-Young procedure)
#' 
#' Inputs are a family of test statistics and a matrix representing a
#' permutation distribution for each member of the family.
#' 
#' The Westfall-Young procedure involves taking the maximum value over family
#' members for each of the N permutations, and comparing all test statistics in
#' the family to this distribution of maximal values.
#' 
#' For recent review in neuroimaging context: https://doi.org/10.1016/j.neuroimage.2020.116760
#' See also: https://dx.doi.org/10.4310/SII.2013.v6.n1.a8
#' 
#' @param tvec A vector of test statistics that constitute a family.
#' @param smat A matrix of simulated test statistics for each family member
#'   acquired through a permutation procedure. Each column corresponds to a test
#'   statistic, and each row a permutation.
#'   
#' @return Adjusted p-values which provide strong control of the FWER.
p_adjust_WestfallYoung <- function(tvec, smat) {
    stopifnot(length(tvec) == ncol(smat))
    p <- apply(smat, MARGIN = 1, FUN = max)
    m <- length(p)
    return(pmin(colSums(outer(p, tvec, '>')) + 1, m) / m)
}
