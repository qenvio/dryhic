#' Apply a bias correction to a HiC contact matrix
#'
#' This function takes a contact matrix and a vector of biases and returns the corrected contact matrix
#' @import Matrix
#' @param mat HiC contact map matrix
#' @param b Vecor of biases
#' @return A corrected HiC contact matrix
#' @details Given a raw contact matrix \eqn{R_{ij}} and a vector of biases \eqn{b_i}, the corrected matrix is \eqn{C_{ij} = R_{ij} / (b_i b_j)}
#' @export
#' @examples
#' plot(0)

correct_mat_from_b <- function(mat, b){

    d <- Diagonal(x = 1 / b)
    out <- d %*% (mat %*% d)
    as(out, class(mat))

    }

