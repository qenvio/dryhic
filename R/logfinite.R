#' Log-transform a vector and substitute infinite values by \code{NA}
#'
#' This function takes a vector of values and returns the log10 transformed values replacing infinite values by \code{NA}
#' @param x A vector of values
#' @return A vector with the log10 transformed values where the infinite ones have been replaced by \code{NA}
#' @export
#' @examples
#' x <- c(10, 1000, 4, 0, 0.1)
#' logfinite(x)

logfinite <- function(x){

    x <- log10(x)
    x[is.infinite(x)] <- NA
    x
    
}
