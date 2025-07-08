library(dplyr)
library(MGnifyR)

## read the data downloaded and stored (see Mgnify.R)
Mgnf <- readRDS("Mgnify_download.RDS")

mg <- MgnifyClient()


Mgnf$studyID <- sapply(Mgnf$relationships$studies$data, function(x){
  paste0(x[["id"]], collapse = "_")
})


Study_tab <- table(Mgnf$studyID)

head(Study_tab[order(Study_tab, decreasing = TRUE)], n=756)
## 756 studies have over 100 samples

sum(head(Study_tab[order(Study_tab, decreasing = TRUE)], n=756))
## These are 361,408 samples!!

## study size will not be a good selection criterion, on the other hand this
## shows that looking through ~700 studies (for metadata etc) will be enough
## to compile a dataset with good sample (health) data


## Post-hoc filtering (as many biome-errors are present)


Mgnf$human_host_tax <- ifelse(Mgnf$attributes$host == "9606",
                              "human", "OTHER")

table(Mgnf$human_host_tax, useNA = "ifany")

Mgnf$biome <- as.character(Mgnf$attributes$`environment-biome`)

table(Mgnf$human_host_tax, Mgnf$attributes$species)
## "\n\nHomo sapiens", "\nHomo sapiens", Home sapiens"
## "Homo sapien" and the correct "Homo sapiens" are wanted

Mgnf$species_human <- ifelse(as.character(Mgnf$attributes$species)%in%
                               c("\n\nHomo sapiens", "\nHomo sapiens",
                                 "Home sapiens", "Homo sapien", "Homo sapiens"),
                             "Hsapiens", "OTHER")

table(Spc = Mgnf$species_human, Host = Mgnf$human_host_tax, useNA = "ifany")


Mgnf$strict_human <- Mgnf$human_host_tax %in% "human" &
  Mgnf$species_human%in%"Hsapiens"


analyses_metadata <- getMetadata(mg, Mgnf$id[[1]])

study_samples <- sapply(unique(Mgnf$studyID)[1:10], function (x) {
  doQuery(mg, "studies", x)
})


Mgnf$intesti <- unlist(lapply(Mgnf$attributes$`sample-metadata`, function (x){
  any(grepl("stool|gut|intesti|fecal|faecal|faeces|feces",
            x[["value"]]))
}))



Mgnf$intesti <- unlist(lapply(Mgnf$attributes$`sample-metadata`, function (x){
  any(grepl("stool|gut|intesti|fecal|faecal|faeces|feces",
            x[["value"]]))
}))

table(Inte = Mgnf$intesti, Hum = Mgnf$strict_human)

### looking a little through metadata
Mgnf[!duplicated(Mgnf$studyID),]$attributes$`sample-metadata`


Mgnf %>%
  group_by(studyID) %>%
  summarize(sample_count = n(),
            human_count = sum(strict_human, na.rm = TRUE),
            envir = paste(unique(biome), collapse = "_")) %>%
  filter(nchar(studyID)>2) %>% arrange(desc(sample_count)) %>%
  print(n=80)


Mgnf$attributes$`sample-metadata`


## export to sort into "good", "bad" and "maybe"
write.csv2(tbl, "biomes_4_manual_annotation.csv")


head(Mgnf)

head(Study_tab[order(Study_tab, decreasing = TRUE)], n=25)

large_studies <- head(Study_tab[order(Study_tab, decreasing = TRUE)], n=25)

large_studies_all <- Mgnf[Mgnf$studyID%in%names(large_studies), ]


large_studies_all$relationships$studies$data

unique(large_studies_all$relationships$studies$data)

tapply(Mgnf$attributes$`environment-biome`, Mgnf$studyID, function(x) table(x))

env_biom_tab <- table(Mgnf$attributes$`environment-biome`)
head(env_biom_tab[order(env_biom_tab, decreasing = TRUE)])

table(Mgnf$attributes$species, useNA = "ifany")

tail(Mgnf$attributes$`sample-metadata`)

table(is.na(Mgnf$attributes$`collection-date`), useNA = "ifany")

loc_tab <- table(Mgnf$attributes$`geo-loc-name`)
head(loc_tab[order(loc_tab, decreasing = TRUE)], n = 50)

table(Mgnf[Mgnf$attributes$`geo-loc-name`%in%"Germany",]$studyID)

Mgnf[Mgnf$attributes$`geo-loc-name`%in%"Germany",]$attributes$`sample-metadata`
## NOT GOOD data!!!

table(Mgnf$attributes$`geo-loc-name`%in%"Malawi")

head(Mgnf$relationships$studies$data)







