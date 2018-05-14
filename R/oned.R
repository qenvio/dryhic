#' Compute oned correction
#'
#' This function takes a \code{data.frame} with bin information and returns a vector of  biases to correct for
#' @importFrom mgcv gam
#' @importFrom mgcv negbin
#' @importFrom mgcv nb
#' @param dat A \code{data.frame} with one bin per row containing the total number of contacts and the potential biases as columns.
#' @param form A \code{formula} describing the  total number of contacts on LHS and the smoothed biases on the RHS (should be compatible with \code{\link{mgcv::gam}})
#' @param p_fit The proportion of data used to fit the model. \code{NA} or \code{NULL} will use the entire data set
#' @param seed Seed number for the sub-sampling (if needed)
#' @param model Logical indicating if model should be returned as \code{attribute} "model"
#' @return A vector of length \code{nrow(dat)} with the biases to correct for.
#' @details Please note that the biases returned are the squarerooted so one can directly apply \code{\link{correct_mat_from_b}}
#' @export
#' @examples
#' plot(0)

oned <- function(dat, form = tot ~ s(map) + s(cg) + s(res),
                 p_fit = NA, seed = 1, model = FALSE){

    if(is.na(p_fit) | is.null(p_fit)){
        dat0 <- dat
    }else{
        set.seed(seed)
        dat0 <- dat[sample(nrow(dat), round(nrow(dat) * p_fit)),]
    }
    
    fit <- mgcv::gam(as.formula(form), data = dat0, family = mgcv::nb())

    out <- predict(fit, newdata = dat, type = "response")

    i_na <- is.na(dat[, all.vars(form)[1]])
    
    out[which(i_na)] <- NA

    out <- as.numeric(sqrt(out / mean(out, na.rm = T)))
    
    if(model){

        attr(out, "model") <- fit

    }

    out

}

