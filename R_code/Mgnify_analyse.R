
Mgnf <- readRDS("Mgnify_download_to5122.RDS")


Mgnf$studyID <- sapply(Mgnf$relationships$studies$data, function(x){
  paste0(x[["id"]], collapse = "_")
})


Study_tab <- table(Mgnf$studyID)

head(Study_tab[order(Study_tab, decreasing = TRUE)], n=254)
## these studies have over 100 samples

sum(head(Study_tab[order(Study_tab, decreasing = TRUE)], n=254))
## Even in a third of the data, this has over 100k samples

## study size will not be a good selection criterion, on the other hand this
## shows that looking through 500-700 studies (for metadata etc) will be enough
## to compile a dataset with good sample (health) data


## Post-hoc filtering (as many biome-errors are present)


Mgnf$human_host_tax <- ifelse(Mgnf$attributes$host == "9606",
                              "human", "OTHER")

table(Mgnf$human_host_tax, useNA = "ifany")

table(Mgnf$human_host_tax, Mgnf$attributes$species)

tbl <- table(Mgnf$attributes$`environment-biome`, Mgnf$human_host_tax)

## export to sort into "good", "bad" and "maybe"
write.csv2(tbl, "biomes_4_manual_annotation.csv")




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







