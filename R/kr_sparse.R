#' Perform the Knight-Ruiz normalization.
#'
#' This function takes a contact matrix and returns its KR normalized version.
#' @param A A \code{dgTMatrix} of contacts with genomic bins as \code{row.names}.
#' @return A A \code{dgTMatrix} of normalized contacts.
#' @details This is a re-implementation of \url{https://github.com/dozmorovlab/HiCcompare/blob/master/R/KRnormalization.R} using sparse matrices.
#' @references \url{https://doi.org/10.1093/imanum/drs019}
#' @export
#' @examples
#' plot(0)

kr_sparse <- function(A) {

    # remove bad rows

    Aold <- A

    b <- attr(A, "b")

    if(!is.null(b)){

        bad <- is.na(b) %>% which
        A <- A[-bad, -bad]

    }

    # initialize

    tol <- 1e-6
    delta <- 0.1
    Delta <- 3
    fl <- 0

    # change NAs in matrix to 0's

    NA_i <- is.na(A@x) %>% which
    A@x[is.na(A@x)] <- 0

    n <- nrow(A)
    e <- matrix(1, nrow = n, ncol = 1) %>% as(class(A))
    x0 <- e
    res <- matrix(nrow = n, ncol = 1) %>% as(class(A))

    # inner stopping criteria

    g <- 0.9
    etamax <- 0.1
    eta <- etamax
    stop_tol <- tol * .5
    x <- x0
    rt <- tol^2
    v <- x *( A %*% x)
    rk <- 1 - v
    rho_km1 <- as.numeric(t(rk) %*% rk)
    rout <- rho_km1
    rold <- rout

    
    MVP <- 0 # We'll count matrix vector products.
    i <- 0  # Outer iteration count.

    while(rout > rt) { # Outer iteration

        i <- i + 1
        k <- 0
        y <- e
        
        innertol <- max(c(eta^2 %*% rout, rt))

        while( rho_km1 > innertol ) { #Inner iteration by CG

            k <- k + 1
            
            if( k == 1) {

                z <- rk / v
                p <- z
                rho_km1 <- as.numeric(t(rk) %*% z)

            }else {

                beta = rho_km1 %*% solve(rho_km2)
                p = z + beta*p
                
            }

            
            # update search direction efficiently

            w <- x * (A %*% (x * p)) + v * p
            alpha <- as.numeric(rho_km1 %*% solve(t(p) %*% w))
            ap <- c(alpha) * p

            # test distance to boundary of cone

            ynew <- y + ap

            if(min(ynew) <= delta) {

                if(delta == 0) break()
                ind <- which(ap < 0)
                gamma <- min((delta - y[ind]) / ap[ind])
                y <- y + gamma %*% ap
                
            }
            
            if(max(ynew) >= Delta) {

                ind <- which(ynew > Delta)
                gamma <- min((Delta-y[ind])/ap[ind])
                y <- y + gamma %*% ap
                
            }
            
            y <- ynew
            rk <- rk - c(alpha) * w
            rho_km2 <- rho_km1
            Z <- rk / v
            rho_km1 <- as.numeric(t(rk) %*% z)
            
        }
        
        x <- x * y
        v <- x*(A %*% x)
        rk <- 1 - v
        rho_km1 <- as.numeric(t(rk) %*% rk)
        rout <- rho_km1
        MVP <- MVP + k + 1

        # Update inner iteration stopping criterion.

        rat <- rout %*% solve(rold)
        rold <- rout
        res_norm <- sqrt(rout)
        eta_o <- eta
        eta <- as.numeric(g %*% rat)

        if(g %*% eta_o^2 > 0.1) {

            eta <- max(c(eta, g %*% eta_o^2))

        }
        
        eta <- max(c(min(c(eta, etamax)), stop_tol / res_norm))
    }

    x <- as.numeric(x)
    
    if(!is.null(b)){

        xx <- b
        xx[!is.na(b)] <- x
        x <- xx
        
    }
    
    result <- Diagonal(x = x) %*% Aold %*% Diagonal(x = x)

    attr(result, "b") <- 1 / x

    result
    
}
