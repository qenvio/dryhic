#' Add pseudo counts to a vector and Log-transform it
#'
#' This function takes a vector of values and returns the log10 transformed values adding a pseudo count to the values prior to the trasformation
#' @param x A vector of values
#' @return A vector with the log10 transformation of the original values plus a pseudo count equal to half of the minimum of the positive original values
#' @details The pseudo count equals half the minimum value, exluding zeros.
#' @export
#' @examples
#' x <- c(10, 1000, 4, 0, 0.1)
#' logpseudo(x)

logpseudo <- function(x) {

    log10(x + min(x[x>0], na.rm = TRUE)/2)

    }
