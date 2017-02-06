#' Compute oned correction
#'
#' This function takes a \code{data.frame} with bin information and returns a vector of  biases to correct for
#' @importFrom mgcv gam
#' @importFrom mgcv negbin
#' @importFrom mgcv nb
#' @param dat A \code{data.frame} with one bin per row containing the total number of contacts and the potential biases as columns.
#' @param form A \code{formula} describing the  total number of contacts on LHS and the smoothed biases on the RHS (should be compatible with \code{\link{mgcv::gam}})
#' @return A vector of length \code{nrow(dat)} with the biases to correct for.
#' @details Please note that the biases returned are the squarerooted so one can directly apply \code{\link{correct_mat_from_b}}
#' @export
#' @examples
#' plot(0)

oned <- function(dat, form = tot ~ s(map) + s(cg) + s(res)){

    fit <- gam(as.formula(form), data = dat, family = nb())

    sqrt(fit$fitted.values / mean(fit$fitted.values))

}

