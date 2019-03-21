#' Compute distance decay of a contact matrix
#'
#' This function takes a contact matrix (tipically the output of \code{\link{get_contacts_matrix}} and returns the distance decay
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr mutate
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
        
        b <- attr(mat, "b")
        
        i_bad <- is.na(b) %>% which
        i_good <- seq(1, nrow(mat)) %>% setdiff(i_bad)

        n_bad <- length(i_bad)
        n_good <- length(i_good)
        n <- n_good + n_bad
        
        mat <- mat[i_good, i_good]

        dec <- mat2df(mat, both = TRUE) %>%
            mutate(d = abs(i_d - j_d),
                   db = abs(i - j),
                   xc = x / (n_good - db)) %>%
            group_by(d) %>%
            summarize(e = sum(xc, na.rm = T)) %>%
            ungroup()
        
        binpos <- rownames(mat) %>%
            gsub("^.*:", "", .) %>%
            as.numeric

        dat <- summary(mat) %>%
            mutate(di = binpos[i],
                   dj = binpos[j],
                   d = abs(di - dj),
                   db = abs(i - j),
                   xc = x / (n_good - db)) %>%
            group_by(d) %>%
            summarize(e = sum(xc, na.rm = T)) %>%
            ungroup()
   
        return(dat)

    }

}
