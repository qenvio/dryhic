#' Compute expected matrix (by distance decay).
#'
#' This function takes a sparse matrix of contacts and data frame with the expected counts per distance and returns a expected matrix.
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @param mat A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}. Currently only one chromosome.
#' @param dd A \code{data.frame} with two columns: \code{d} (distance in bp) and \code{e} (expected number of contacts.
#' @return A \code{dgTMatrix} corresponding to expected counts.
#' @export
#' @examples
#' plot(0)

get_expected_matrix <- function(mat, dd){

    if(!(class(mat) == "dgTMatrix" | class(mat) == "dgCMatrix")) stop("mat should be a sparse matrix")

    bad_bins <- attr(mat, "b") %>% is.na() %>% which()

    out <- mat2df(mat, both = TRUE) %>%
        mutate(d = abs(i_d - j_d)) %>%
        left_join(dd, by = "d") %>%
        with(sparseMatrix(i = i, j = j, x = x / e,
                          dims = dim(mat),
                          dimnames = dimnames(mat))) %>%
        as(class(mat))

    attr(out, "b") <- attr(mat, "b")
    
    out

}
