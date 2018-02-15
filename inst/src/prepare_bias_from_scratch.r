#!/usr/bin/env Rscript   

# dependencies

library("dplyr")

set_remote_file <- function(genome){
    paste0("https://bismap.hoffmanlab.org/raw/",
           genome,
           "/k50.umap.bedgraph.gz")
}

options(stringsAsFactors = F)

# settings

genome <- "hg19"
motif <- "GATC"
resolution <- as.integer(1e5)
n_cores <- 2
species <- "Hsapiens"

# load genome package

library(paste0("BSgenome.", species, ".UCSC.", genome), character.only = T)
genome_ucsc <- get(species)

chromosomes <- names(genome_ucsc) %>% grep("_|chrM", ., value = T, invert = T)

# get chromosome sizes and bins

chromosome_sizes <- sapply(chromosomes, function(x) genome_ucsc[[x]] %>% length) %>%
    data.frame(chr = names(.),
               size = .)

bins <- group_by(chromosome_sizes, chr) %>%
    do(data.frame(chr = .$chr,
                  pos = seq(0,
                            floor(.$size / resolution) * resolution, by = resolution) %>%
                      as.integer)) %>%
    ungroup

# get RE coordinates

re_coordinates <- parallel::mclapply(chromosomes,
                          function(x) {
                              out <- matchPattern(motif, genome_ucsc[[x]], fixed = F)
                              idx <- as.character(out)
                              idx <- idx == motif
                              out <- out[idx,]
                              data.frame(chr = x, start = start(out), end = end(out))
                          }, mc.cores = n_cores, mc.preschedule = F) %>%
    do.call(rbind, .)

re_bins <- mutate(re_coordinates,
                  pos = as.integer(floor(start / resolution) * resolution)) %>%
    group_by(chr, pos) %>%
    summarize(re = n()) %>%
    ungroup %>%
    right_join(bins, by = c("chr", "pos"))

re_bins$re[is.na(re_bins$re)] <- 0

# get CG content

cg_bins <- lapply(chromosomes, function(chrom){
    chrseq <- genome_ucsc[[chrom]]
    idx <- seq(1, length(chrseq) - resolution, by = resolution)
    nuc_freq <- letterFrequencyInSlidingView(chrseq,
                                             resolution,
                                             c("A", "C", "G", "T", "N"))[idx,]
    last_bin <- seq(floor(length(chrseq) / resolution) * resolution,
                    length(chrseq))
    last_bin <- letterFrequency(chrseq[last_bin],
                                c("A", "C", "G", "T", "N"))
    nuc_freq <- rbind(nuc_freq, last_bin)
    rownames(nuc_freq) <- NULL
    nuc_freq <- filter(bins, chr == chrom) %>%
        cbind(nuc_freq)
    mutate(nuc_freq,
           cg = (C + G) / (A + C + G + T + N)) %>%
        select(chr, pos, cg)
}) %>%
    do.call(rbind, .)


# get mappability

tmp <- tempfile()

set_remote_file(genome) %>%
    download.file(tmp)

map_bed <- read.delim(tmp, skip = 1, head = F)
names(map_bed) <- c("chr", "start", "end", "map")


map_bins <- mutate(map_bed, size = end - start,
                   pos = as.integer(floor(start / resolution) * resolution)) %>%
    group_by(chr, pos) %>%
    summarize(map = weighted.mean(map, size)) %>%
    ungroup %>%
    right_join(bins, by = c("chr", "pos"))

map_bins$map[is.na(map_bins$map)] <- 0


# merge all info

bias <- inner_join(re_bins,
                   cg_bins,
                   by = c("chr", "pos")) %>%
    inner_join(map_bins,
               by = c("chr", "pos"))


