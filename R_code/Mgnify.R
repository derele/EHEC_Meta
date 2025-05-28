library(httr)
library(jsonlite)
library(dplyr)
library(purrr)


base_url <- "https://www.ebi.ac.uk/metagenomics/api/v1/samples"
biome <- "root:Host-associated:Human:Digestive system"
url <- paste0(base_url, "?biome=", URLencode(biome, reserved = TRUE))

samples <- NULL
count <- 0

repeat {
  res <- GET(url)
  j <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
  new_samples <- as_tibble(j$data)

  if (is.null(samples)) {
    samples <- new_samples
  } else {
    samples <- bind_rows(samples, new_samples)
  }

  count <- count + 1
  url <- j$links[["next"]]
  message("Finished round: ", count)

  if (url==j$links$last) break
}

#### on the Server (because of Problems on some pages):

## Finished round: 5122

##  Fehler: lexical error: invalid char in json text.
##                                       <!doctype html> <html lang="en"
##                     (right here) ------^


## saveRDS(samples, "Mgnify_download_to5122.RDS")

##  > url <- j$links[["next"]]
##  > url
## [1] "https://www.ebi.ac.uk/metagenomics/api/v1/samples?biome=root%3AHost-associated%3AHuman%3ADigestive+system&page=5123"


##  restart repeat loop without resetting "count" or overwriting "samples"!

### AGAIN ## AGAIN ## AGAIN

saveRDS(samples, "Mgnify_download_to_7402.RDS")

## url <- j$links[["next"]]
## url

## [1] "https://www.ebi.ac.uk/metagenomics/api/v1/samples?biome=root%3AHost-associated%3AHuman%3ADigestive+system&page=7403"

##  restart repeat loop without resetting "count" or overwriting "samples"!

saveRDS(samples, "Mgnify_download.RDS")




