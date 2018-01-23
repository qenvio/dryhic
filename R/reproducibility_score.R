#' Compute the reproducibility score between two contact matrices.
#'
#' This function takes a couple of contact matrices and computes the reproducibility score.
#' @param mat1 A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param mat2 Another \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @param r The number of componentes used in the comparison.
#' @param mipr Threshold for the j-inverse participation ratio (IPR).
#' @return The reproducibility score between two matrices.
#' @details Please note that mat1 and mat2 are expected to be a matrices of one chromosome (you could obtain them using \code{subset_matrix}).
#' @references \url{https://doi.org/10.1093/bioinformatics/btx152}
#' @export
#' @examples
#' plot(0)

reproducibility_score <- function(mat1, mat2, r = 20, mipr = 5){

    na1 <- summary(mat1) %>%
        filter(is.na(x)) %>%
        (function(x) c(x$i, x$j) %>% unique)

    na2 <- summary(mat2) %>%
        filter(is.na(x)) %>%
        (function(x) c(x$i, x$j) %>% unique)

    if(length(na1) > 0){
        
        mat1 <- mat1[-na1, -na1]

    }

    if(length(na2) > 0){

        mat2 <- mat2[-na2, -na2]

    }

    if(min(nrow(mat1), nrow(mat2)) <= 1) return(NA)
    
    bins <- intersect(rownames(mat1), rownames(mat2))
    b1 <- match(bins, rownames(mat1))
    b2 <- match(bins, rownames(mat2))

    mat1 <- mat1[b1, b1]
    mat2 <- mat2[b2, b2]
    
    i1 <- symmetrize_matrix(mat1) %>% rowSums(na.rm = T)
    i2 <- symmetrize_matrix(mat2) %>% rowSums(na.rm = T)

    i1 <- which(i1 > 0)
    i2 <- which(i2 > 0)

    i <- intersect(i1, i2)

    mat1 <- mat1[i,i]
    mat2 <- mat2[i,i]
    
    l1 <- laplacian_matrix(mat1)
    l2 <- laplacian_matrix(mat2)
    
    r <- min(r, nrow(l1)) %>% min(nrow(l2))

    ei1 <- RSpectra::eigs(l1, r, which = "SM")
    ei2 <- RSpectra::eigs(l2, r, which = "SM")
   
    ev1 <- matrix(0, nrow(l1), r)
    ev1 <- Re(ei1$vectors[, 1:r])
    ev2 <- matrix(0, nrow(l2), r)
    ev2 <- Re(ei2$vectors[, 1:r])

    ipr1 <- 1 / colSums(ev1 ^ 4)
    ipr2 <- 1 / colSums(ev2 ^ 4)

    ev1 <- ev1[, ipr1 > mipr, drop = F]
    ev2 <- ev2[, ipr2 > mipr, drop = F]

    if((ncol(ev1) == 0) | (ncol(ev2) == 0)) return(NA)
    
    r <- min(ncol(ev1), ncol(ev2))

    out_sd <- sapply(1:r, function(i) evec_dist(ev1[,i], ev2[,i])) %>% sum

    abs(sqrt(2) - out_sd / r)/sqrt(2)

}
