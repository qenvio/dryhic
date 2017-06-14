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
    
    oldreso <- rownames(mat) %>%
        strsplit(":") %>%
        do.call(rbind, .) %>%
        as.data.frame %>%
        setNames(c("chr", "pos")) %>%
        mutate(pos = as.integer(pos)) %>%
        group_by(chr) %>%
        arrange(pos) %>%
        summarize(m = diff(pos) %>% min) %>%
        ungroup  %>%
        summarize(m = min(m)) %$%
        m

    if(oldreso == newreso) return(mat)

    if((oldreso %% newreso) != 0) stop("New resolution should be multiple of the original one")
    
    ids <- rownames(mat)
    chrs <- gsub(":.*$", "", ids)
    pos <- gsub("^.*:", "", ids) %>% as.integer

    newpos <- as.integer(floor(pos / newreso) * newreso)

    newids <- paste(chrs, newpos, sep = ":")

    m <- fac2sparse(newids)

    newmat <- m %*% (mat %*% t(m))

    bins <- rownames(mat)[rownames(mat) %in% rownames(newmat)]

    i <- match(bins, rownames(newmat))
    
    as(newmat[i,i], class(mat))

}
