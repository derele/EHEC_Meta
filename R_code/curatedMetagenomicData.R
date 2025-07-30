library(curatedMetagenomicData)
library(dplyr)
library(magrittr)

## fully anotated disease status!!!
table(is.na(sampleMetadata$disease))

table(sampleMetadata$body_site)
## great stool      21030

sampleMetadata %>%
  filter(body_site == "stool") %>%
  pull(NCBI_accession) %>%
  writeLines("data/cMD_ncbi.txt")
