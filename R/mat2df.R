#' Transform the contacts of a \code{dgTMatrix} into a \code{data.frame}
#'
#' This function takes a sparse matrix of contacts and transforms it into a \code{data.frame}. It works only with matrices of one chromosome.
#' @param mat A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param b An optional \code{vector} of genomic position (otherwhise, it uses the information contained in the \code{row.names}.
#' @param both A \code{logical} indicating if both (bin index and position) should be returned.
#' @return A \code{data.frame} with the two coordinates (\code{i} and \code{j}) and the number of contacts (\code{x}).
#' @details Please note that mat is expected to be a matrix of one chromosome (you could obtain one using \code{subset_matrix}).
#' @export
#' @examples
#' plot(0)

mat2df <- function(mat, b = NULL, both = FALSE){

    if(!(class(mat) == "dgTMatrix" | class(mat) == "dgCMatrix")) stop("mat should be a sparse matrix")
    
    if(is.null(b)) {

        b <- rownames(mat) %>%
            gsub("^.*:", "", .) %>%
            as.integer

    }

    if(both) {

        out <- summary(mat) %>%
            mutate(i_d = b[i],
                   j_d = b[j])

    }else {
        
        out <- summary(mat) %>%
            mutate(i = b[i],
                   j = b[j])

    }

    out

}
