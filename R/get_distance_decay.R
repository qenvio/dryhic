#' Compute distance decay of a contact matrix
#'
#' This function takes a contact matrix (tipically the output of \code{\link{get_contacts_matrix}} and returns the distance decay
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @param mat A contact matrix (could be the output of \code{\link{get_contacts_matrix}})
#' @param reso Bin size (resolution) in bp
#' @return If \code{mat} contains more than one chromosome, a list of \code{data.frame}s. If not, a \code{data.frame} with the distance \code{d} (in bp) and the expected number of counts \code{e}.
#' 
#' @export
#' @examples
#' plot(0)

get_distance_decay <- function(mat, reso){

    chroms <- rownames(mat) %>% gsub(":.*$", "", .) %>% unique()

    if(length(chroms) > 1){

        dat <- lapply(chroms, function(chrom){
            subset_matrix(mat, chrom) %>%
                get_distance_decay(reso)
        }) %>%
            setNames(chroms)

        return(dat)
        
    }else{
        
        n <- nrow(mat)
        dec <- mat2df(mat, both = T) %>%
            filter(i >= j) %>%
            mutate(d = i_d - j_d,
                   dd = i - j) %>%
            group_by(d, dd) %>%
            summarize(s = sum(x, na.rm = T),
                      na = is.na(x) %>% sum()) %>%
            ungroup() %>%
            mutate(e = s / (n - dd - na)) %>%
            dplyr::select(d, e)
        return(dec)

    }

}
