library(curatedMetagenomicData)
library(dplyr)
library(magrittr)
library(tidyverse)
library(pheatmap)


write_NEW_SRA <- FALSE
write_2nd_SRA <- FALSE

## fully anotated disease status!!!
table(is.na(sampleMetadata$disease))

table(sampleMetadata$body_site)
## great stool      21030

sM <- sampleMetadata %>%
  filter(body_site == "stool")


if(write_NEW_SRA){
  sM %>%
    filter(!is.na(NCBI_accession)) %>%
    pull(NCBI_accession) %>%
    writeLines("data/cMD_ncbi.txt")
}

accs <- sM %>%
  pull(NCBI_accession)


## Extracting multiple entries in NCBI Accession column
multiAcc <- unique(accs[!grepl("^\\w{3}\\d{6,7}$", accs)])
multiAcc_list <- sapply(multiAcc, strsplit, ";")


multiAcc_second <- lapply(multiAcc_list, function(x) x[2:length(x)])
missingAcc <- unlist(multiAcc_second)
missingAcc <- missingAcc[!is.na(missingAcc)]

if(write_2nd_SRA){
    writeLines(missingAcc, "data/cMD_2nd_ncbi.txt")
}

### reading the data after processing with "Other_code/map_sra_bowDia.sh"
results_dir <- "H:/Analysenetz/transfer/EHEC_processed/"

bt_raw <- list.files(results_dir, pattern = "\\.bt2\\.idxstats\\.txt$", full.names = TRUE) %>%
  map(~ {
    tab <- read.table(.x, sep = "\t", stringsAsFactors = FALSE)
    tab$sample <- sub("\\.bt2\\.idxstats\\.txt$", "", basename(.x))
    tab
  }) %>%
  bind_rows() %>%
  as_tibble() %>%
  rename(
    gene_id = V1,
    length = V2,
    mapped = V3,
    unmapped = V4
  )

bt_collated <- bt_raw %>%
  filter(gene_id != "*") %>%
  mutate(marker = sub("(.*)\\(gb\\|.*;(.*)\\)", "bt_\\1_\\2", gene_id)) %>%
  select(sample, marker, mapped) %>%
  pivot_wider(
    names_from = marker,
    values_from = mapped,
    values_fill = 0
  )

table(STX2A = bt_collated$bt_VFG000837_stx2A, STX2B = bt_collated$bt_VFG000838_stx2B)
table(bt_collated$bt_VFG000739_eae>1, STX2B = bt_collated$bt_VFG000803_eae>1)


dia_files <- list.files(results_dir, pattern="*.diamond.tsv", full.names = TRUE)

summarise_diamond <- function(file) {
  sample_id <- str_remove(basename(file), "\\.diamond\\.tsv$")

  if (file.size(file) == 0) {
    return(tibble(
      marker        = "VFG000739(gb|AAC38392;eae)", ## any one as placeholder
      hits          = 0L,
      hits_95  = 0L,
      sample        = sample_id
    ))
  }

  df <- read_tsv(file,
                 col_names = c(
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"
  ), show_col_types = FALSE)

  df %>%
    group_by(marker = sseqid) %>%
    summarise(
      hits = n(),
      hits_95 = sum(pident > 95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_id)
}

## Apply to all files and combine
dia_summary <- bind_rows(lapply(dia_files, summarise_diamond))

dia_collated <- dia_summary %>%
  mutate(marker = sub("(.*)\\(gb\\|.*;(.*)\\)", "\\1_\\2", marker)) %>%
  select(sample, marker, hits, hits_95) %>%
  pivot_wider(
    names_from = marker,
    values_from = c(hits, hits_95),
    values_fill = 0
  )


pheatmap(log10(dia_collated[, !colnames(dia_collated)%in%"sample"]+1))
dia_collated95 <- dia_collated[, grepl("hits_95", colnames(dia_collated))]
pheatmap(log10(dia_collated95[rowSums(dia_collated95)>2, ]+1))


Screen <- inner_join(bt_collated, dia_collated)
Screen_df <- as.data.frame(Screen)

rownames(Screen_df) <- Screen_df$sample

specificRS <- rowSums(Screen[, grepl("bt_|hits_95_", colnames(Screen))])

stx2RS <- rowSums(Screen[, grepl("stx2", colnames(Screen))])

pheatmap(log10(Screen_df[stx2RS>1, grepl("stx2", colnames(Screen_df))]+1),
        )

Stx2_samples <- rownames(Screen_df[stx2RS>1,])

controls <- c("SRR12195337", "SRR12195338", "SRR12195339",
              "SRR4101313", "SRR4101314", "SRR4101315", "SRR33438544",
              "SRR32314265")

length(Stx2_samples)/nrow(Screen)*100

sM[sM$NCBI_accession%in%Stx2_samples&
  !sM$NCBI_accession%in%controls   , 1:25]


### need to expand the NCBI_ID here once mutliple runs per samples are read!

finalData <- inner_join(sM, Screen, by = c("NCBI_accession" = "sample"))

foo <- finalData[finalData$NCBI_accession%in%Stx2_samples, ]


