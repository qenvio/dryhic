#' Bin the genome
#'
#' This function takes a HiC-BAM file and creates a data.frame with the genome binned at the desired resolution
#' @import magrittr
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @param inbam HiC-BAM file (only the header will be used)
#' @param resolution Desired resolution (bin size)
#' @return A \code{data.frame} containing chromosome, position and bin ID
#' @export
#' @examples
#' plot(0)

make_bins <- function(inbam, resolution){

    sizes <- paste("samtools view -H", inbam) %>%
        pipe %>%
        read.delim(head = F) %>%
        filter(V1 == "@SQ") %>%
        mutate(V2 = gsub("^SN:", "", V2),
               V3 = gsub("^LN:", "", V3)) %>%
        (function(x) as.numeric(x$V3) %>% setNames(x$V2))

    bins <- lapply(sizes, function(x) seq(0, x, resolution) %>% as.integer)

    bins <- data.frame(chr = rep(names(bins), sapply(bins, length)),
                       pos = unlist(bins),
                       stringsAsFactors = F) %>%
        mutate(bin = paste(chr, pos, sep = ":"))

    bins
    
}
