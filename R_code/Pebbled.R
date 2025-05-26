library(SRAdb)
library(ggplot2)
library(dplyr)
library(xml2)
library(rentrez)

sqlfile <- "/data/db/SRAmetadb.sqlite"
sra_con <- dbConnect(RSQLite::SQLite(), sqlfile)

peb <- read.delim("pebblescout-meta-summary.tsv")

length(peb$SubjectID)
length(unique(peb$SubjectID))

## how many (unique) Samples per Study (proxy: Title)
ttab <- tapply(peb$BioSample, peb$Title, function(x) length(unique(x)))
ttab <- ttab[order(ttab, decreasing = TRUE)]

head(ttab)

## Now more robust and correct
## How many (unique) Samples per Study (per Study ID)
## How many were part of the Study, how many positive

biosample_ids <- unique(peb$BioSample)

## Starting to work with SRA

dbListTables(sra_con)  # List all tables
dbListFields(sra_con, "study")  # Check columns in the study table
dbListFields(sra_con, "sample")  # Check columns in the sample table
dbListFields(sra_con, "metaInfo")

formatted_ids <- paste0("'", paste(biosample_ids, collapse = "','"), "'")

filtered_query <- sprintf("
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


filtered_result <- dbGetQuery(sra_con, filtered_query)

filtered_result$short_title <- substr(filtered_result$study_title, 1, 40)
filtered_result$prevalence_raw <- filtered_result$matching_biosample_count/
    filtered_result$total_biosample_count

## filtered_result[ , !colnames(filtered_result)%in%c("study_title",
##                                                    "study_abstract")]

filtered_result$twin_study <- grepl("twin", filtered_result$study_abstract, ignore.case=TRUE)

filtered_result$healthy_study <- grepl("healthy", filtered_result$study_abstract,
                                       ignore.case=TRUE)

filtered_result <- filtered_result %>%
    mutate(species_study = case_when(
               grepl("patient|adult|human|twin|child|infant|participant",
                     study_abstract, ignore.case = TRUE) ~ "human",
               grepl("bovine|rumen|calve", study_abstract, ignore.case = TRUE) ~ "bovine",
               grepl("sewage|manure", study_abstract, ignore.case = TRUE) ~ "environmenal",
               TRUE ~ "unknown"  # Default case
           ))


table(filtered_result$study_type, filtered_result$species_study)


filtered_result_with_labels <- filtered_result[filtered_result$total_biosample_count > 250, ]

ggplot(filtered_result, aes(y = prevalence_raw, x = species_study,
                            color = twin_study,
                            shape = study_type,
                            size = total_biosample_count)) +
    geom_jitter() ## +
    ## geom_text(data = filtered_result_with_labels,
    ##           aes(label = short_title),
    ##           vjust = -0.5,  # Adjust vertical position
    ##           hjust = 1,     # Adjust horizontal position
    ##           check_overlap = TRUE)  # Avoid overlapping text



### We need to do better stuff on the sample-level!

## taxon_id "408170" human gut metagenome is our TARGET
## taxon_id "1504969" human blood metagenome is BAD


## Many samples are missing information or have incomplete information
## we also seem to lose samples

detailed_query <- sprintf("
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


detailed_query <- sprintf("
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


detailed_result <- dbGetQuery(sra_con, detailed_query)
head(detailed_result)
nrow(detailed_result)

table(biosample_ids%in%detailed_result$sample_alias, useNA="ifany")

missingIDs <- biosample_ids[!biosample_ids%in%detailed_result$sample_alias]

f_missing_ids <- paste0("'", paste(missingIDs, collapse = "','"), "'")

missing_query <- sprintf("
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
", f_missing_ids)


## still missing!
## head(missing_result)


missing_query <- sprintf("
  SELECT
    sample.sample_ID,
    sample.sample_alias,
    sample.sample_accession,
    study.study_ID,
    study.study_accession
  FROM sample
  LEFT JOIN study ON sample.submission_accession = study.submission_accession
  WHERE sample.sample_alias IN (%s) OR sample.sample_accession IN (%s)
", f_missing_ids, f_missing_ids)

missing_result <- dbGetQuery(sra_con, missing_query)

## still missing
head(missing_result)


dbListTables(sra_con)  # List all tables
dbListFields(sra_con, "study")  # Check columns in the study table
dbListFields(sra_con, "sample")  # Check columns in the sample table
dbListFields(sra_con, "experiment")  # Check columns in the experiment table



Entrez_list <- lapply(missingIDs[1:8], function(ID) {
    result <- entrez_search(db = "biosample", term = ID)
    entrez_summary(db = "biosample", id = result$ids)
})


parse_biosample <- function(sample_data) {
    parsed_xml <- read_xml(sample_data)
    tibble(
        biosample_id = xml_text(xml_find_first(parsed_xml, "//Id[@db='BioSample']")),
        sample_name = xml_text(xml_find_first(parsed_xml, "//Id[@db_label='Sample name']")),
        sra_id = xml_text(xml_find_first(parsed_xml, "//Id[@db='SRA']")),
        title = xml_text(xml_find_first(parsed_xml, "//Description/Title")),
        taxonomy_id = xml_attr(xml_find_first(parsed_xml, "//Organism"), "taxonomy_id"),
        taxonomy_name = xml_attr(xml_find_first(parsed_xml, "//Organism"), "taxonomy_name"),
        owner = xml_text(xml_find_first(parsed_xml, "//Owner/Name")),
        first_name = xml_text(xml_find_first(parsed_xml, "//Contacts/Contact/Name/First")),
        last_name = xml_text(xml_find_first(parsed_xml, "//Contacts/Contact/Name/Last")),
        model = xml_text(xml_find_first(parsed_xml, "//Models/Model")),
        package = xml_text(xml_find_first(parsed_xml, "//Package")),
        host = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='host']")),
        isolation_source = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='isolation_source']")),
        collection_date = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='collection_date']")),
        geo_loc_name = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='geo_loc_name']")),
        lat_lon = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='lat_lon']")),
        sample_processing = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='samp_mat_process']")),
        sample_size = xml_text(xml_find_first(parsed_xml, "//Attribute[@attribute_name='samp_size']")),
        bioproject_id = xml_text(xml_find_first(parsed_xml, "//Link[@type='entrez' and @target='bioproject']"))
    )
}

## Apply the function to all samples
all_biosamples <- lapply(Entrez_list, parse_biosample)

## Combine the individual tibbles into a single table
biosample_table <- bind_rows(all_biosamples)

head(biosample_table)
