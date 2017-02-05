#' Plot HiC contact matrix
#'
#' This function takes a HiC contact matrix and creates heatmap-like prepresentation of it
#' @import magrittr
#' @import Matrix
#' @param mat HiC contact matrix (it could be the output of \code{\link{get_contacts_matrix}})
#' @param coord A vector of size two with the start and end coordinates of the desired region to plot
#' @param resolution Resolution (bin size) in bp
#' @param transformation Transformation to apply to the contacts prior to the heatmap representation
#' @param color Vector of colors
#' @param sym Logical indicating if the color scale should be symmetrical (around 0)
#' @param trim Numeric proportion of the data that should be 'flattened' on both extremes
#' @param rotate Logical indicating if the matrix should be rotated
#' @param unit_x_axis Numeric unit for the X axis
#' @param label_x_axis Character X axis label
#' @param na.col Color to depict \code{NA} values
#' @param ... Further arguments to be passed to \code{\link{plot}}
#' @export
#' @examples
#' plot(0)

plot_matrix <- function(mat, coord, resolution,
                        transformation = logfinite,
                        color = colorRampPalette(c("white", "red"))(100),
                        sym = FALSE, trim = .01, rotate = FALSE,
                        unit_x_axis = 1e6,
                        label_x_axis = "Genomic Position / Mbp",
                        na.col = "white",
                        ...){

    # settings

    options(scipen = 999)
    
    # prepare matrix

    rownames(mat) <- colnames(mat) <- rownames(mat) %>% gsub("^.*:", "", .)
    
    lims <- seq(floor(coord[1] / resolution),
                ceiling(coord[2] / resolution)) * resolution
    
    i <- rownames(mat) %in% lims
    
    mat <- as.matrix(mat[i, i])
    mat <- mat[match(lims, rownames(mat)), match(lims, rownames(mat))]
    rownames(mat) <- colnames(mat) <- lims

    if(rotate) mat <- mat[nrow(mat):1,]
    
    # prepare axis info and parameters

    guides <- pretty(x = rownames(mat) %>% as.numeric)

    guides_pos <- data.frame(y = 1:nrow(mat), x = rownames(mat) %>% as.numeric) %>%
        lm(y ~ x, .) %>%
        predict(newdata = data.frame(x = guides))
    
    par(mar = c(4, 0, 0, 0), pty = "s")

    # trimm

    if(trim > 0){

        trim <- as.matrix(mat) %>% c %>% quantile(c(trim / 2, 1 - trim / 2), na.rm = T)
        mat[mat < trim[1]] <- trim[1]
        mat[mat > trim[2]] <- trim[2]

    }
    
    # prepare range of colors

    x <- as.matrix(mat) %>% transformation %>% as.matrix
    if(sym){

        upper <- max(abs(x), na.rm = T)
        lower <- - upper
        
    }else{
        
        lower <- min(x, na.rm = T)
        upper <- max(x, na.rm = T)
        
    }

    # transform scores into colors

    if(max(x, na.rm = T) == min(x, na.rm = T)){
        x[] <- color[round(length(color) / 2)]
    }else{
        x[] <- color[cut(c(x), seq(lower, upper,
                                   len = length(color) + 1), include = T)]
    }

    x[is.na(x)] <- na.col
    
    # get limits of genomic region

    range_pos <- as.numeric(rownames(x)) %>% range

    # get matrix dimensions
    
    nr <- nrow(x)
    nc <- ncol(x)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d

    # plot void region

    plot(NA, type="n",
         xlim=c(0, nr), ylim=c(0, nc),
         xlab=label_x_axis,
         ylab="",
         asp=1, axes = F, cex.lab = 1.5, ...)

    # add heatmap

    rasterImage(as.raster(unclass(x)),
                xleft = 0, xright = nc, ybottom = 0, ytop = nr,
                interpolate = FALSE)
    axis(1, at = guides_pos,
         labels = guides / unit_x_axis, cex.axis = 1.5)
    
    invisible()
    
}
