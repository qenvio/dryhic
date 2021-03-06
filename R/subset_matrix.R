#' Extract the contacts associated to a region form a acontact matrix.
#'
#' This function takes a sparse matrix of contacts and a genomic region and extract a sub-matrix with the desired region.
#' @param mat A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param region A genomic region ("chr" or "chr:start-end").
#' @return A \code{dgTMatrix} corresponding to the desired region.
#' @export
#' @examples
#' plot(0)

subset_matrix <- function(mat, region){

    if(!(class(mat) == "dgTMatrix" | class(mat) == "dgCMatrix")) stop("mat should be a sparse matrix")

    chrom <- gsub(":.*$", "", region)

    bins <- rownames(mat)

    bins <- data.frame(chr = gsub(":.*$", "", bins),
                       pos = gsub("^.*:", "", bins) %>% as.numeric,
                       i = 1:length(bins))

    if(!grepl("-", region)){
        i <- filter(bins, chr == chrom) %$% i
    }else{
        rr <- gsub("^.*:", "", region) %>% strsplit("-") %>% unlist() %>% as.numeric()
        i <- filter(bins, chr == chrom, pos >= rr[1], pos <= rr[2]) %$% i
    }

    out <- mat[i,i]

    b <- attr(mat, "b")

    if(!is.null(b)){

        attr(out, "b") <- b[i]

    }

    out

}
