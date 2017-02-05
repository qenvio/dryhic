#' Get upper triangular matrix from a HiC contact matrix
#'
#' This function takes a HiC contact matrix and returns the upper triangular matrix. It divides the diagonal counts by 2.
#' @import magrittr
#' @import Matrix
#' @param mat HiC contact matrix
#' @return The upper triangular matrix of the original one
#' @export
#' @examples
#' plot(0)

upper_tri <- function(mat){

    check <- c(triu(mat, -1) %>% sum(na.rm = T), tril(mat, -1) %>% sum(na.rm = T))
    
    if(min(check) == 0) {
        
        message("Matrix seems to be already upper diagonal, doing nothing ...")
        
        return(mat)

    }

    out <- triu(mat)

    as(out - band(out, 0, 0) / 2, class(mat))

    }
