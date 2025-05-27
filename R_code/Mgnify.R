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

saveRDS(samples, "Mgnify_download.RDS")

samples$studyID <- sapply(samples$relationships$studies$data, function (x) {
  paste(x[["id"]], collapse = "_")
})

table(samples$studyID)

### Wow! This data is so so rich!
