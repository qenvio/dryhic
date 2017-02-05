#' Quick and dirty implementation of ICE matrix normalization
#'
#' This function takes a contact matrix (tipically the output of \code{\link{get_contacts_matrix}} and performs ICE normalization
#' @import magrittr
#' @import Matrix
#' @importFrom dplyr mutate
#' @param mat A contact matrix (could be the output of \code{\link{get_contacts_matrix}})
#' @param intermax Maximum number of iterations
#' @param maxdev Maximum relative deviation from the mean of the correction factor 
#' @param verbose Logical indicating if progression messages should be turned on
#' @return A named \code{sparse_matrix} containing the normalized number of contacts per pair of genomic bins at the requested region plus the \code{b} attribute containing the vector of correction factors
#' @references @article{imakaev2012iterative,
#'   title={Iterative correction of Hi-C data reveals hallmarks of chromosome organization},
#'   author={Imakaev, Maxim and Fudenberg, Geoffrey and McCord, Rachel Patton and Naumova, Natalia and Goloborodko, Anton and Lajoie, Bryan R and Dekker, Job and Mirny, Leonid A},
#'   journal={Nature methods},
#'   volume={9},
#'   number={10},
#'   pages={999--1003},
#'   year={2012},
#'   publisher={Nature Publishing Group}
#' }
#' 
#' @export
#' @examples
#' plot(0)

ICE <- function(mat, itermax = 1, maxdev = .1, verbose = TRUE){

    b <- attr(mat, "b")
    
    if(is.null(b)){
           
        b <- rep(1, ncol(mat))

    }

    for(i in 1:itermax){
    
        s <- colSums(mat, na.rm = T)
        s[is.na(b)] <- NA
        sm <- mean(s, na.rm = T)
        
        delta_b <- (s/sm)
        
        B <- Diagonal(x = 1/delta_b)
        
        b <- b * delta_b
        
        mat <- B %*% (mat %*% B)
        
        dev <- abs((s/sm - 1)) %>% max(na.rm = T)
        
        if (verbose) {
            
            paste(i,
                  min(s, na.rm = T) %>% format(digits = 2),
                  sm %>% format(digits = 2),
                  max(s, na.rm = T) %>% format(digits = 2),
                  dev %>% format(digits = 5),
                  "\n",
                  sep = "\t") %>%
                cat
        }
        
        if(dev < maxdev) break

    }

    attr(mat, "b") <- b
    
    mat

}
