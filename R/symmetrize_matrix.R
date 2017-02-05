#' Get the full symmetric matrix from a upper triangular contact matrix
#'
#' This function takes a upper triangular HiC contact matrix and returns the full symmetric matrix.
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr summarize
#' @importFrom dplyr slice
#' @param mat HiC contact matrix (should be upper triangular)
#' @return The full symmetric HiC contact matrix corresponding to the input one
#' @export
#' @examples
#' plot(0)

symmetrize_matrix <- function(mat){

    check <- c(triu(mat, -1) %>% sum(na.rm = T), tril(mat, -1) %>% sum(na.rm = T))

    if(min(check) > 0) {
        
        message("Matrix seems to be already symmetric, doing nothing ...")
        
        return(mat)

    }

    mat <- mat + t(mat)

    mat

}
