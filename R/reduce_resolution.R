#' Reduce resolution (incrrease bin size) of a HiC contact matrix
#'
#' This function takes a contact matrix and returns the corresponding contact matrix with the resolution reducction (increase of bin size)
#' @import magrittr
#' @import Matrix
#' @param mat HiC contact map matrix
#' @param newreso Desired bin size in bp
#' @return A HiC contact matrix with the desired resolution
#' @export
#' @examples
#' plot(0)

reduce_resolution <- function(mat, newreso){

    ids <- rownames(mat)
    chrs <- gsub(":.*$", "", ids)
    pos <- gsub("^.*:", "", ids) %>% as.numeric

    newpos <- floor(pos / newreso) * newreso

    newids <- paste(chrs, newpos, sep = ":")

    m <- fac2sparse(newids)

    newmat <- m %*% (mat %*% t(m))

    colnames(newmat) <- rownames(newmat) <- rownames(newids)

    as(newmat, class(mat))

}
