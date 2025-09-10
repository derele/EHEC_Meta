library(curatedMetagenomicData)
library(dplyr)
library(magrittr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(boot)


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


map_results <- inner_join(bt_collated, dia_collated)
map_results_df <- as.data.frame(map_results)

rownames(map_results_df) <- map_results_df$sample

## need to work on the controls!!
controls <- c("SRR12195337", "SRR12195338", "SRR12195339",
              "SRR4101313", "SRR4101314", "SRR4101315", "SRR33438544",
              "SRR32314265")


### Some samples have the same NCBI accession as they are the same used in
## multiple studies (shoud have downloaded and screened them only once; TODO:
## select unique accessions when completing the dataset!!)

### we here discard simple the second mention of ever accession in the data!
## they come from currated data!

sM <- sM[!duplicated(sM$NCBI_accession), ]

## now we need to expand the multiple collected accessions.

library(dplyr)

# Beispiel: sM_long hat die Spalten: group_id (alte Zeilennummer) + NCBI_accession
sM_long <- sM %>%
  mutate(my_sample_id = row_number()) %>%        # jede alte Zeile bekommt eine ID
  separate_rows(NCBI_accession, sep = ";")   # aufsplitten

# Nun joinen wir mit der Screen-Tabelle
joined <- sM_long %>%
  left_join(map_results, by = c("NCBI_accession" = "sample"))

# Jetzt für jede Gruppe aufsummieren:
collated_joined <- joined %>%
  group_by(my_sample_id) %>%
  summarise(
    across(starts_with("bt_") | starts_with("hits_"), sum, na.rm = TRUE),
    across(!starts_with("bt_") & !starts_with("hits_"), ~ first(.x)),
    .groups = "drop"
  )

## need to expand the NCBI_ID here once mutliple runs per samples are read!
finalData <- inner_join(sM, map_results, by = c("NCBI_accession" = "sample")) %>%
  as.data.frame()

rownames(finalData) <- finalData$my_sample_ID

finalData <- collated_joined
rownames(finalData) <- finalData$my_sample_id

screenCols <- colnames(finalData[, grepl("bt_|hits_", colnames(finalData))])
specificCols <- colnames(finalData[, grepl("bt_|hits_95_", colnames(finalData))])
stxCols <- colnames(map_results[, grepl("stx", colnames(map_results))])
stx2Cols <- colnames(map_results[, grepl("stx2", colnames(map_results))])

### back to the contreols in the overall screening dataset... to see some
## non-human (non curatedMetagenomcData merging)
pheatmap(log10(map_results_df[controls , screenCols]+1))


posStudies <- rownames(finalData)[rowSums(finalData[, screenCols])>5]
pposStudies <- rownames(finalData)[rowSums(finalData[, specificCols])>5]
stxStudies <- rownames(finalData)[rowSums(finalData[, stxCols])>3]
stx2Studies <- rownames(finalData)[rowSums(finalData[, stx2Cols])>0]

finalData$DE2011_outbreak <- ifelse(finalData$study_name=="LomanNJ_2013",
                                    "Yes", "No")

study_annot <- finalData[, c("DE2011_outbreak", "disease", "age", "study_name")]
study_annot$age <- as.numeric(study_annot$age)

pheatmap(log10(finalData[posStudies, screenCols]+1),
##         annotation_row = study_annot[, c(1, 3)],
         show_rownames = FALSE)

## based on this (and controls we should drop eae < 95%)
screenCols <- screenCols[!screenCols%in%c("hits_VFG000739_eae", "hits_VFG000803_eae")]

pheatmap(log10(finalData[pposStudies, screenCols]+1),
##         annotation_row = study_annot[, c(1, 3)],
         show_rownames = FALSE)

pheatmap(log10(finalData[stxStudies, screenCols]+1),
##         annotation_row = study_annot[, c(1, 2, 3)],
         show_rownames = FALSE)

pheatmap(log10(finalData[stx2Studies, screenCols]+1),
##         annotation_row = study_annot[, c(2,3)],
         show_rownames = FALSE)


pheatmap(log10(finalData[finalData$disease%in%"STEC", screenCols]+1),
        ## annotation_row = study_annot[, c(2,3)],
         show_rownames = FALSE)


length(Stx2_samples)/nrow(finalData)*100

tapply(finalData$number_bases, finalData$study_name, mean)/1000000

table(finalData[!finalData$study_name%in%"LomanNJ_2013", "bt_VFG000837_stx2A"])
table(finalData[!finalData$study_name%in%"LomanNJ_2013",
                "bt_VFG000837_stx2A"]>0)

finalData$STX2bt <- finalData$bt_VFG000837_stx2A>0|
  finalData$bt_VFG000838_stx2B>0

nrow(finalData[finalData$STX2bt & !finalData$study_name%in%"LomanNJ_2013", ])/
  nrow(finalData[!finalData$study_name%in%"LomanNJ_2013", ])*100
### 0.3624883%

## 0.280095%

## but detection is dependent on sequencing depth
depth_model<- glm(STX2bt ~ number_bases,
                 data=finalData[finalData$disease%in%"healthy",],
                 family = binomial)

## From LomanNJ_2013 depth at their 0.67% sensitivity
# 1. Get reference prediction
ref_depth <- data.frame(number_bases = 1.075e9)
current_prob <- predict(depth_model, newdata = ref_depth, type = "response")

target_prob <- 0.67
current_linear_pred <- predict(depth_model, newdata = ref_depth, type = "link")  # Log-odds
beta_1 <- coef(depth_model)["number_bases"]
new_intercept <- log(target_prob / (1 - target_prob)) - beta_1 * ref_depth$number_bases

depth_model_adjusted <- depth_model
depth_model_adjusted$coefficients["(Intercept)"] <- new_intercept

predict(depth_model_adjusted, newdata = ref_depth, type = "response")  # Should output 0.67

new_data <- data.frame(number_bases = seq(1e8, 4e10, length.out = 100))
new_data$prob <- predict(depth_model_adjusted, newdata = new_data, type = "response")



finalData[finalData$disease%in%"healthy", "pred_sensitivity"] <-
  predict(depth_model_adjusted,
          new_data = finalData[finalData$disease%in%"healthy","number_bases"],
          type="response") + 0.67


# 4. Plot
ggplot(new_data, aes(x = number_bases, y = prob)) +
  geom_line() +
  geom_vline(xintercept = 1.075e9, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.67, linetype = "dashed", color = "red") +
  # geom_point(data=finalData, aes(x = number_bases, y = pred_sensitivity,
  #                                color=study_name)) +
  labs(x = "Sequencing Depth (number_bases)", y = "Adjusted Detection Probability",
       title = "Logistic Model Calibrated to 67% Sensitivity at 1.075B Reads")



ggsave("figures/detection_depht.png")



## At about 4e+10 = 40,000,000,000 == 40 GB we can expect to have 100% sensitivity.

table(finalData$number_bases > 4e+10) ## no 100% sensitivity
table(finalData$number_bases > 1.3e+10) ## 188 samples with 95% sensitivity at 13GB data

table(finalData[!finalData$study_name%in%"LomanNJ_2013", "STX2bt"])

finalData$pred_sensitivity <- predict(depth_model_adjusted,
                                      newdata = finalData,
                                      type = "response")

true_positives <- sum(finalData$STX2bt / finalData$pred_sensitivity)

true_positives/nrow(finalData)*100

set.seed(123)
boot_fn <- function(data, indices) {
  sample_data <- data[indices, ]
  obs_pos <- sum(sample_data$STX2bt)
  pred_sens <- predict(depth_model_adjusted, newdata = sample_data, type = "response")
  true_pos <- sum(sample_data$STX2bt / pred_sens)
  true_prev <- true_pos / nrow(sample_data)
  return(true_prev)
}
boot_results <- boot(finalData, boot_fn, R = 1000)
quantile(boot_results$t, c(0.025, 0.975))  # 95% CI

boot_df <- data.frame(true_prevalence = boot_results$t)

ggplot(boot_df, aes(x = true_prevalence)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30,
                 fill = "skyblue",
                 color = "black",
                 alpha = 0.7) +
  geom_density(linewidth = 1, color = "darkblue") +
  geom_vline(xintercept = mean(boot_results$t),
             linetype = "solid",
             color = "red",
             linewidth = 1) +
  geom_vline(xintercept = quantile(boot_results$t, c(0.025, 0.975)),
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  labs(
    x = "Prevalence of stx2AB",
    y = "Density",
    title = "Bootstrap Distribution of Prevalence Estimates",
    subtitle = "Red lines: Mean (solid) and 95% CI (dashed)"
  ) +
  theme_minimal()


ggsave("figures/bootstrap_estimate.png")


disease_model<- glm(STX2bt ~ disease,
                    ##                  data=finalData[!finalData$study_name%in%"LomanNJ_2013",],
                    data=finalData,
                    family = binomial)

summary(disease_model)
