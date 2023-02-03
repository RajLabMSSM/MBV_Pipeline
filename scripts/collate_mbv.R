#Jack Humphrey 2023
# collate MBV outputs and make one nice table with best matches


library(tidyverse)
library(optparse)


option_list <- list(
        make_option(c('--outFolder'), help='', default = "."),
        make_option(c('--outFile'), help='', default = "test.tsv")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder
outFile <- opt$outFile

# get all bamstat outputs
files <- list.files(path = outFolder, pattern = "*bamstat.txt", full.names = TRUE)

# read in and pick best matching genotype for each rna-seq sample
pick_top <- function(data){
    sample <-  gsub(".bamstat.txt", "", basename(data) )
    d <- readr::read_delim(data, delim = " ", col_types = "cnnnnnnnnnn")
    top <- dplyr::arrange(d, desc(perc_het_consistent) ) %>% 
        head(1) %>% 
        dplyr::mutate(sample_id = sample) %>%
        dplyr::select(sample_id, genotype_id = SampleID, everything() ) 
    return(top)
}

all <- purrr::map_df(files, pick_top)

readr::write_tsv(all, outFile)



 
