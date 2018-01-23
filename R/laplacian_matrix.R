#' Get the associated laplacian matrix.
#'
#' This function takes sparse matrix of contacts and returns its corresponding Laplacian matrix.
#' @param x A \code{dgTMatrix} of contacts.
#' @return The associated \code{dgTMatrix} laplacian matrix.
#' @details This is an auxiliary function used by \code{reproducibility_score}.
#' @export
#' @examples
#' plot(0)

laplacian_matrix <- function(x){

    if(!(class(x) == "dgTMatrix" | class(x) == "dgCMatrix")) stop("mat should be a sparse matrix")
    
    x <- symmetrize_matrix(x)
    d <- rowSums(x)

    D <- Diagonal(x = d)

    L <- D - x

    D <- Diagonal(x = 1 / sqrt(d))
    
    (D %*% L) %*% D

}
