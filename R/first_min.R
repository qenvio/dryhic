#' Find first local minimum
#'
#' This funciton computes the lowest minimum value of the density estimate of theinput numeric vector
#' @importFrom dplyr lag
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @param x A numeric \code{vector} of values.
#' @param plot \code{logical} indicating if histogram should be plotted.
#' @param zero_rm \code{logical} indicating if zeros should be removed from calculation.
#' @return A numeric value of the lower local minimum value of the density estimation.
#' @export
#' @examples
#' plot(0)

first_min <- function(x, plot = FALSE, zero_rm = FALSE){

    x <- na.omit(x)

    if(zero_rm){
        x <- x[x > 0 & x < quantile(x, .99)]
    }else{
        x <- x[x < quantile(x, .99)]
    }

    out <- density(x) %>%
        with(data.frame(x = x, y = y)) %>%
        mutate(d = y - lag(y),
               s = sign(d),
               dd = s - lag(s)) %>%
        filter(x > 0, dd == 2) %>%
        arrange(x) %>%
        head(1) %$%
        x %>%
        round()

    if(plot){
        hist(x, 100, las = 1)
        abline(v = out, col = "red")
    }

    out
    
}
