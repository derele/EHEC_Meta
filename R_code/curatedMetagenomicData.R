library(curatedMetagenomicData)
library(dplyr)
library(magrittr)
library(tidyverse)

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
    pull(NCBI_accession) %>%
    writeLines("data/cMD_ncbi.txt")
}

accs <- sM %>%
  pull(NCBI_accession)

multiAcc <- unique(accs[!grepl("^\\w{3}\\d{6,7}$", accs)])
multiAcc_list <- sapply(multiAcc, strsplit, ";")


multiAcc_second <- lapply(multiAcc_list, function(x) x[2:length(x)])

missingAcc <- unlist(multiAcc_second)

if(write_2nd_SRA){
    writeLines(missingAcc, "data/cMD_2nd_ncbi.txt")
}

### reading the data after processing with "Other_code/map_sra_bowDia.sh"
results_dir <- "H:/Analysenetz/transfer/EHEC_processed/"

bar <- list.files(results_dir, pattern = "\\.bt2\\.idxstats\\.txt$", full.names = TRUE) %>%
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

baz <- bar %>%
  filter(gene_id != "*") %>%
  mutate(marker = sub("(.*)\\(gb\\|.*;(.*)\\)", "bt_\\1_\\2", gene_id)) %>%
  select(sample, marker, mapped) %>%
  pivot_wider(
    names_from = marker,
    values_from = mapped,
    values_fill = 0
  )

table(STX2A = baz$bt_VFG000837_stx2A, STX2B = baz$bt_VFG000838_stx2B)
table(baz$bt_VFG000739_eae>1, STX2B = baz$bt_VFG000803_eae>1)


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

dia_baz <- dia_summary %>%
  mutate(marker = sub("(.*)\\(gb\\|.*;(.*)\\)", "\\1_\\2", marker)) %>%
  select(sample, marker, hits, hits_95) %>%
  pivot_wider(
    names_from = marker,
    values_from = c(hits, hits_95),
    values_fill = 0
  )

library(pheatmap)

pheatmap(log10(dia_baz[, !colnames(dia_baz)%in%"sample"]+1))
dia_baz95 <- dia_baz[, grepl("hits_95", colnames(dia_baz))]
pheatmap(log10(dia_baz95[rowSums(dia_baz95)>2, ]+1))


Screen <- inner_join(baz, dia_baz)
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


