library(SRAdb)
library(ggplot2)
library(dplyr)
library(xml2)
library(rentrez)
library(magrittr)
library(dplyr)
library(tidyr)
library(pheatmap)

sqlfile <- "H:/Analysenetz/transfer/SRAmetadb.sqlite"
sra_con <- dbConnect(RSQLite::SQLite(), sqlfile)

## Pebblescout was used setting "Max #subjects per query" to 1,000,000
## "Score constant" at default value of 13,000 and "Masking threshold" at
## default value of 5663


peb_files <- list.files("data/", pattern = "pebblescout-meta-summary_*")

peb <- lapply(peb_files, function(x) {
  cbind(read.delim(paste0("data/", x)),
        gene=gsub("pebblescout-meta-summary_(.*)\\.tsv", "\\1", x))
}) %>% do.call(rbind, .) %>% as_tibble()

## number of hits
length(peb$SubjectID)
## number of genes hit
length(unique(peb$SubjectID))

tapply(peb$SubjectID, peb$gene, function(x) length(unique(x)))

tapply(peb$SubjectID, peb$gene,
       function(x) length(unique(x)))

tapply(peb$BioSample, peb$gene,
       function(x) length(unique(x)))

all.genes <- as.factor(unique(peb$gene))

peb <- peb[nchar(peb$BioSample) > 0, ]


peb[peb$gene%in%"VirDb", "gene"] <- substr(as.vector(peb[peb$gene%in%"VirDb", "QueryID"])[[1]], 1, 4)

COpebbled <- peb %>%
  mutate(presence = TRUE) %>%
  distinct(BioSample, gene, presence) %>%
  pivot_wider(
    names_from = gene,
    values_from = presence,
    values_fill = FALSE
  )



pheatmap(apply(COpebbled[, -1], 2, as.numeric))

string_peb <- peb[peb$PBSscore>95, ]

sCOpebbled <- string_peb %>%
  mutate(presence = TRUE) %>%
  distinct(BioSample, gene, presence) %>%
  pivot_wider(
    names_from = gene,
    values_from = presence,
    values_fill = FALSE
  )

pheatmap(apply(sCOpebbled[, -1], 2, as.numeric))


## how many (unique) Samples per Study (proxy: Title)
ttab <- tapply(peb$BioSample, peb$Title, function(x) length(unique(x)))
ttab <- ttab[order(ttab, decreasing = TRUE)]

head(ttab)

## Now more robust and correct
## How many (unique) Samples per Study (per Study ID)
## How many were part of the Study, how many positive

biosample_ids <- unique(peb$BioSample)

sbiosample_ids <- unique(string_peb$BioSample)


## Starting to work with SRA

dbListTables(sra_con)  # List all tables
dbListFields(sra_con, "study")  # Check columns in the study table
dbListFields(sra_con, "sample")  # Check columns in the sample table
dbListFields(sra_con, "metaInfo")

formatted_ids <- paste0("'", paste(biosample_ids, collapse = "','"), "'")

sformatted_ids <- paste0("'", paste(sbiosample_ids, collapse = "','"), "'")


study_query <- sprintf("
  SELECT
    study.study_accession AS bioproject,
    study.study_title AS study_title,
    study.study_type AS study_type,
    study.study_abstract as study_abstract,
    study.study_description as study_description,
    study.study_attribute as study_attribute,
    COUNT(DISTINCT CASE WHEN sample.sample_alias IN (%s) THEN sample.sample_alias END) AS matching_biosample_count,
    COUNT(DISTINCT sample.sample_alias) AS total_biosample_count
  FROM sample
  JOIN study ON sample.submission_accession = study.submission_accession
  GROUP BY study.study_accession, study.study_title
  HAVING matching_biosample_count > 0
", formatted_ids)

sstudy_query <- sprintf("
  SELECT
    study.study_accession AS bioproject,
    study.study_title AS study_title,
    study.study_type AS study_type,
    study.study_abstract as study_abstract,
    study.study_description as study_description,
    study.study_attribute as study_attribute,
    COUNT(DISTINCT CASE WHEN sample.sample_alias IN (%s) THEN sample.sample_alias END) AS matching_biosample_count,
    COUNT(DISTINCT sample.sample_alias) AS total_biosample_count
  FROM sample
  JOIN study ON sample.submission_accession = study.submission_accession
  GROUP BY study.study_accession, study.study_title
  HAVING matching_biosample_count > 0
", sformatted_ids)


study_result <- dbGetQuery(sra_con, study_query)

sstudy_result <- dbGetQuery(sra_con, sstudy_query)


## We lose many
sum(study_result$matching_biosample_count)
length(biosample_ids)
length(unique(peb$BioSample))

study_result$short_title <- substr(study_result$study_title, 1, 40)
study_result$prevalence_raw <- study_result$matching_biosample_count/
    study_result$total_biosample_count

study_result$twin_study <- grepl("twin", study_result$study_abstract, ignore.case=TRUE)

study_result$healthy_study <- grepl("healthy", study_result$study_abstract,
                                       ignore.case=TRUE)

study_result <- study_result %>%
    mutate(species_study = case_when(
               grepl("patient|adult|human|twin|child|infant|participant",
                     study_abstract, ignore.case = TRUE) ~ "human",
               grepl("bovine|rumen|calve", study_abstract, ignore.case = TRUE) ~ "bovine",
               grepl("sewage|manure", study_abstract, ignore.case = TRUE) ~ "environmenal",
               TRUE ~ "unknown"  # Default case
           ))

sstudy_result$prevalence_raw <- sstudy_result$matching_biosample_count/
  sstudy_result$total_biosample_count




### We need to do better stuff on the sample-level!

## taxon_id "408170" human gut metagenome is our TARGET
## taxon_id "1504969" human blood metagenome is BAD


## Many samples are missing information or have incomplete information
## we also seem to lose samples

sample_query <- sprintf("
  SELECT
    sample.sample_ID,
    sample.sample_alias,
    sample.sample_accession,
    study.study_ID AS study_ID,
    study.study_accession AS study_acc,
    study.study_alias AS study_ali
  FROM sample
  JOIN study ON sample.submission_accession = study.submission_accession
  WHERE sample.sample_alias IN (%s)
", formatted_ids)


sample_query <- sprintf("
  SELECT
    sample.sample_ID,
    sample.sample_alias,
    sample.sample_accession,
    sample.taxon_id,
    study.study_ID AS study_ID,
    study.study_accession AS study_acc,
    study.study_alias AS study_ali,
    study.submission_accession AS study_sub_acc,
    sample.submission_accession AS sample_sub_acc
  FROM sample
  LEFT OUTER JOIN study
    ON sample.submission_accession = study.submission_accession
  WHERE sample.sample_alias IN (%s)
", formatted_ids)


sample_result <- dbGetQuery(sra_con, sample_query)
head(sample_result)
nrow(sample_result)

table(biosample_ids%in%sample_result$sample_alias, useNA="ifany")

sum(study_result$matching_biosample_count)

table(sample_result$study_acc%in%study_result$bioproject, useNA="ifany")
table(study_result$bioproject%in%sample_result$study_acc, useNA="ifany")
length(unique(sample_result$study_acc))
### ... okay this fits 165 studies!

## have a manual look!!
study_result[order(study_result$prevalence_raw, decreasing=TRUE),
             c("total_biosample_count", "short_title")]


write.csv2(study_result, file="manual_pebble.csv")

study_result <- read.csv2(file="manual_pebble_FINISHED.csv")


write.csv2(study_result, file="manual_pebble.csv")



table(study_result$study_type, study_result$verdict)


study_result %>%
  filter(grepl("Good", study_result$verdict)) %>%
  ggplot(aes(x=total_biosample_count, y = prevalence_raw,
                         size = total_biosample_count,
             color = verdict)) +
  geom_point() +
  scale_y_continuous("raw prevalence %", breaks = c(1:4/10),
                     labels = c(1:4*10)) +
  scale_x_continuous("total numbe of samples in study")



### below I try to find the missing IDs without success... commented in the
## interest of running this here faster.

# missingIDs <- biosample_ids[!biosample_ids%in%sample_result$sample_alias]
#
# f_missing_ids <- paste0("'", paste(missingIDs, collapse = "','"), "'")
#
# missing_query <- sprintf("
#   SELECT
#     sample.sample_ID,
#     sample.sample_alias,
#     sample.sample_accession,
#     sample.taxon_id,
#     study.study_ID AS study_ID,
#     study.study_accession AS study_acc,
#     study.study_alias AS study_ali,
#     study.submission_accession AS study_sub_acc,
#     sample.submission_accession AS sample_sub_acc
#   FROM sample
#   LEFT OUTER JOIN study
#     ON sample.submission_accession = study.submission_accession
#   WHERE sample.sample_alias IN (%s)
# ", f_missing_ids)

#
# missing_result <- dbGetQuery(sra_con, missing_query)
#
#
# ## still missing
# head(missing_result)

# missing_query <- sprintf("
#   SELECT
#     sample.sample_ID,
#     sample.sample_alias,
#     sample.sample_accession,
#     study.study_ID,
#     study.study_accession
#   FROM sample
#   LEFT JOIN study ON sample.submission_accession = study.submission_accession
#   WHERE sample.sample_alias IN (%s) OR sample.sample_accession IN (%s)
# ", f_missing_ids, f_missing_ids)
#
# missing_result <- dbGetQuery(sra_con, missing_query)
#
# ## still missing
# head(missing_result)


