#' Compute the stratum adjusted correlation coefficient (SCC) between two contact matrices.
#'
#' This function takes a couple of contact matrices and computes the reproducibility score.
#' @param mat1 A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param mat2 Another \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param mm The maximum distance (in bp) to truncate the data.
#' @param resol The matrix resolution (bin size).
#' @return A list with the estimated SCC and associated std.
#' @details Please note that mat1 and mat2 are expected to be a matrices of one chromosome (you could obtain them using \code{subset_matrix}).
#' @references \url{https://doi.org/10.1101/101386 }
#' @export
#' @examples
#' plot(0)

scc <- function(mat1, mat2, mm, resol){

    unitize <- function(x){
        rank(x, ties.method = "random") / length(x)
    }

    dat <- paste_matrices(mat1, mat2, mm) %>% na.omit
    
    if(min(unique(dat$x1) %>% length,
           unique(dat$x2) %>% length) == 1) {
        
        return(list(scc = NA, std = NA))
        
    }
    
    res <- mutate(dat, d = abs(i - j) / resol) %>%
        dplyr::group_by(d) %>%
        mutate(r1 = unitize(x1),
               r2 = unitize(x2)) %>%
        summarize(corr = cor(x1, x2),
                  wei = sqrt(var(r2) * var(r1)) * n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(-d) %>%
        na.omit %>%
        as.list

    res$scc <- with(res, corr * wei/sum(wei)) %>% as.numeric %>% sum
    res$std <- with(res, sqrt(sum(wei^2 * var(corr))/(sum(wei))^2)) %>% as.numeric

    res

}
