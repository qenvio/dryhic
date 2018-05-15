#' Create a \code{data.frame} with the contacts of two matrices (including zeros).
#'
#' This function takes two sparse matrix of contacts and a maximum distance (in bins) and returns a \code{data.frame} with the contacts of the tow matrices for each pair of positions.
#' @param smat1 A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param smat2 Another \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param mm The maximum distance (in bp) to truncate the data.
#' @return A \code{data.frame} with the two coordinates (\code{i} and \code{j}) and the number of contacts of each matrix (\code{x1} and \code{x2}).
#' @details This is an auxiliary function used by \code{scc}.
#' @export
#' @examples
#' plot(0)

paste_matrices <- function(smat1, smat2, mm){

    if(!(class(smat1) == "dgTMatrix" | class(smat1) == "dgCMatrix")) stop("smat1 should be a sparse matrix")

        if(!(class(smat2) == "dgTMatrix" | class(smat2) == "dgCMatrix")) stop("smat2 should be a sparse matrix")

    out1 <- mat2df(smat1) %>%
        as.data.frame

    p1 <- rownames(smat1) %>% gsub("^.*:", "", .) %>% as.numeric

    out1 <- mutate(out1, x1 = x, x = NULL)

    out2 <- mat2df(smat2) %>%
        as.data.frame

    p2 <- rownames(smat2) %>% gsub("^.*:", "", .) %>% as.numeric

    out2 <- mutate(out2, x2 = x, x = NULL)

    pb <- intersect(p1, p2)

    scaf <- expand.grid(i = pb, j = pb) %>%
        filter(j - i <= mm, i <= j)
    
    out <- dplyr::full_join(out1, out2, by = c("i", "j")) %>%
        dplyr::right_join(scaf, by = c("i", "j"))

    out[is.na(out)] <- 0
    
    out
    
}
