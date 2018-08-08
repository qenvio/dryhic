#' Compute the stratum adjusted correlation coefficient (SCC) between two contact matrices.
#'
#' This function takes a couple of contact matrices and computes the reproducibility score.
#' @import magrittr
#' @import data.table
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
    unitize <- function(x) {
        rank(x, ties.method = "random")/length(x)
    }
    if (!(class(smat1) == "dgTMatrix" | class(smat1) == "dgCMatrix")) 
        stop("smat1 should be a sparse matrix")
    if (!(class(smat2) == "dgTMatrix" | class(smat2) == "dgCMatrix")) 
        stop("smat2 should be a sparse matrix")

    out1 <- mat2df(smat1) %>% as.data.frame
    p1 <- rownames(smat1) %>% gsub("^.*:", "", .) %>% as.numeric
    out1 <- mutate(out1, x1 = x, x = NULL) %>%
        filter(j - i <= mm, i <= j)

    out2 <- mat2df(smat2) %>% as.data.frame
    p2 <- rownames(smat2) %>% gsub("^.*:", "", .) %>% as.numeric
    out2 <- mutate(out2, x2 = x, x = NULL) %>%
        filter(j - i <= mm, i <= j)

    pb <- intersect(p1, p2)

    aux <- bandSparse(nrow(smat1), ncol(smat1),
                      k = 0:(mm / resol),
                      diagonals = matrix(1,
                                         nrow(smat1),
                                         mm / resol + 1)) %>%
        as(class(smat1))
    rownames(aux) <- colnames(aux) <- rownames(smat1)

    scaf <- mat2df(aux) %>%
        filter(i %in% pb, j %in% pb) %>%
        dplyr::select(-x)

    out1 <- data.table(out1)
    out2 <- data.table(out2)
    scaf <- data.table(scaf)

    setkey(out1, i, j)
    setkey(out2, i, j)
    setkey(scaf, i, j)

    out <- out1[out2[scaf]]

    out[is.na(x1), x1 := 0]
    out[is.na(x2), x2 := 0]

    out[, d:= abs(i - j)/resol]
    
    if (min(unique(out$x1) %>% length, unique(out$x2) %>% length) == 
        1) {
        return(list(scc = NA, std = NA))
    }

    out <- out[, .(r1 = unitize(x1),
                   r2 = unitize(x2),
                   x1,
                   x2),
               .(d)]
    out <- out[, .(corr = cor(x1, x2),
                   wei = sqrt(var(r2) * var(r1)) * .N),
               .(d)]

    res <- as.list(out)
    res$scc <- with(res,
                    corr * wei/sum(wei)) %>%
        as.numeric() %>% 
        sum()
    res$std <- with(res, (sum(wei^2 * var(corr))/(sum(wei))^2) %>%
                         sqrt()) %>% 
        as.numeric()

    res[c("corr", "wei", "scc", "std")]

}
