#' Make a contact matrix from HiC data
#'
#' This function takes a HiC-BAM file and creates a contact map as a \code{sparse_ matrix} with \code{dimnames} corresponding to the genomic bins
#' @import Matrix
#' @import magrittr
#' @importFrom dplyr mutate
#' @param inbam HiC-BAM file
#' @param resolution Desired resolution (bin size)
#' @param pos A vector with the genomic position IDs to be included in the output. It can be the \code{bin} column of the output from \code{\link{make_bins}}
#' @param region Region of interest following the format \code{chr|chr:start-end}. If \code{NULL}, take the entire genome
#' @param whole Logical indicating if the entire available genome should be retrieved (rectangular contact map) or only the region of interest indicated with \code{region} (square contact map). It should be consistent the \code{pos} argument
#' @param filtin Integer value indicating the required bits for the HiC-BAM FLAG (include the corresponding contacts in the output)
#' @param filtex Integer value indicating the excluded bits for the HiC-BAM FLAG (exclude the corresponding contacts in the output)
#' @return A named \code{sparse_matrix} containing the number of contacts per pair of genomic bins at the requested region
#' @export
#' @examples
#' plot(0)

get_contacts_matrix <- function(inbam, resolution, pos, region = NULL,
                                whole = F, filtin = 0, filtex = 783){

    if(whole){
        
        script_file <- system.file("src", "bam_to_mat_whole.sh", package = "dryhic")

    }else{

        script_file <- system.file("src", "bam_to_mat.sh", package = "dryhic")
        
    }
    
    paste(script_file,
          "-f", filtin,
          "-F", filtex,
          "-w", resolution,
          inbam,
          region) %>%
        pipe %>%
        read.delim(head = F, stringsAsFactors = F) %>%
        mutate(b1 = factor(paste(V1, V2, sep = ":"), levels = pos),
               b2 = factor(paste(V3, V4, sep = ":"), levels = pos)) %>%
        xtabs(V5 ~ b1 + b2, ., sparse = T)
    
}
