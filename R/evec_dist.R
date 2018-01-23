#' Compute the distance between two eigenvectors.
#'
#' This function computes the euclidian distance between two eigenvectros, taking into account their arbitrary relative signs.
#' @param ev1 A eigen vector.
#' @param ev2 Another eigen vector.
#' @return The euclidean distance between the vectors..
#' @export
#' @examples
#' plot(0)

evec_dist <- function(ev1, ev2){

    d1 <- sum((ev1 - ev2)^2) %>% sqrt
    d2 <- sum((ev1 + ev2)^2) %>% sqrt

    ifelse(d1 < d2, d1, d2)

}
