library(readr)
library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(purrr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

fname <- args[1]
outfile <- args[2]

d <- read_tsv(fname, 
    col_types = cols_only(
        anc_base = 'c',
        mut_base = 'c',
        mut_type = 'c')) %>%
  filter(mut_type == 'SNP') %>%
  select(contains('base')) %>%
  group_by(anc_base, mut_base) %>%
  tally()

write_tsv(d, outfile)
