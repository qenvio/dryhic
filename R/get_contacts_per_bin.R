#' Retrieve HiC data per genomic bin
#'
#' This function takes a HiC-BAM file and creates a nammed \code{vector} with the total number of contact per genomic bin
#' @import magrittr
#' @importFrom dplyr mutate
#' @param inbam HiC-BAM file
#' @param resolution Desired resolution (bin size)
#' @param region Region of interest following the format \code{chr|chr:start-end}. If \code{NULL}, take the entire genome
#' @param pos A vector with the genomic position IDs to be included in the output. It can be the \code{bin} column of the output from \code{\link{make_bins}}
#' @param filtin Integer value indicating the required bits for the HiC-BAM FLAG (include the corresponding contacts in the output)
#' @param filtex Integer value indicating the excluded bits for the HiC-BAM FLAG (exclude the corresponding contacts in the output)
#' @return A named \code{vector} containing the total number of contacts per genomic bin at the requested region
#' @export
#' @examples
#' plot(0)

get_contacts_per_bin <- function(inbam, resolution, region = NULL,
                                 pos, filtin = 0, filtex = 783){
    
    script_file <- system.file("src", "bam_to_bin.sh", package = "dryhic")
    
    out <- paste(script_file,
                 "-f", filtin,
                 "-F", filtex,
                 "-w", resolution,
                 inbam,
                 region) %>%
        pipe %>%
        read.delim(head = F, stringsAsFactors = F) %>%
        mutate(b1 = factor(paste(V1, V2, sep = ":"), levels = pos)) %>%
        xtabs(V3 ~ b1, .)
    
}
