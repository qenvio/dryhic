#' Subsample a HiC contact matrix
#'
#' This function takes a HiC contact matrix and returns a subsample of it. The input matrix should be triangular.
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr summarize
#' @importFrom dplyr slice
#' @param mat HiC contact map matrix
#' @param p Numeric probability of being sub-sampled
#' @return A sub-sampled HiC contact matrix
#' @export
#' @examples
#' plot(0)

subsample_matrix <- function(mat, p = .5){

    check <- c(triu(mat, -1) %>% sum(na.rm = T), tril(mat, -1) %>% sum(na.rm = T))

    if(min(check) != 0) stop("Matrix should be triangular")
    
    mat@x <- rbinom(length(mat@x), mat@x, prob = p) %>% as.numeric

    drop0(mat)

}
