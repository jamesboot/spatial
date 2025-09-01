# Script to annotate CARD results with spatial meta data

setwd('Documents/Milner_et_al/')

# Packages
library(tidyverse)

# Load spatial data
partekMeta <- read.delim('exported.txt')

# Locate CARD outputs
cardFiles <- list.files(path = 'CARD_outs',
                        pattern = '*.csv',
                        full.names = T)

# Create sample names from file locs
samples <- gsub('CARD_outs/', '', gsub('.csv', '', cardFiles))

# Load CARD outputs and name list elements
cardDatList <- lapply(cardFiles, function(x) {
  read.csv(x)
})
names(cardDatList) <- samples

# Change col names to allow easy merging later
cardDatList <- lapply(cardDatList, function(x) {
  colnames(x)[1] <- 'Cell.name'
  x
})

# Now merge for all samples in for loop
for (x in samples) {
  # Subset down to sample of interest
  subset <- partekMeta %>%
    filter(Sample.name == x)
  # Merge
  merged <-
    list(cardDatList[[x]], subset) %>% purrr::reduce(left_join, by = 'Cell.name')
  # Write to file
  write.csv(merged, file = paste0('CARD_outs/', x, '_plusMeta.csv'))
}
