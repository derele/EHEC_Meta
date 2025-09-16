library(curatedMetagenomicData)
library(dplyr)
library(magrittr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(boot)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

write_NEW_SRA <- FALSE
write_2nd_SRA <- FALSE
write_final_SRA <-FALSE

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

map_results <- full_join(bt_collated, dia_collated)
## to figure out which samples are missing just one of the methods!
map_results_all <- inner_join(bt_collated, dia_collated)

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


# my_sample_id (alte Zeilennummer) + NCBI_accession
sM_long <- sM %>%
  mutate(my_sample_id = row_number()) %>%        # jede alte Zeile bekommt eine ID
  separate_rows(NCBI_accession, sep = ";")   # aufsplitten

## we use map_results_all because we want to complete the collection using data
### from both mapping methods!!
table(sM_long$NCBI_accession%in%map_results_all$sample)

stillMissingACC <- sM_long$NCBI_accession[!sM_long$NCBI_accession%in%map_results_all$sample]

if(write_final_SRA){
  writeLines(stillMissingACC, "data/cMD_3rd_ncbi.txt")
}

# Nun joinen wir mit der Screen-Tabelle
joined <- sM_long %>%
  left_join(map_results, by = c("NCBI_accession" = "sample"))

# Jetzt für jede Gruppe aufsummieren:
finalData <- joined %>%
  group_by(my_sample_id) %>%
  summarise(
    across(starts_with("bt_") | starts_with("hits_"), sum, na.rm = TRUE),
    across(!starts_with("bt_") & !starts_with("hits_"), ~ first(.x)),
    .groups = "drop"
  )

## bowtie but super sensitive setting everything to stx2 when detected
finalData$STX2bt <- finalData$bt_VFG000837_stx2A>0|
  finalData$bt_VFG000838_stx2B>0

## this is the positive control: finalData$study_name=="LomanNJ_2013"
finalData$DE2011 <- ifelse(finalData$study_name=="LomanNJ_2013",
                           "YES", "NO")

## some samples have zero sequencing depth
finalData <- subset(finalData, number_bases > 0)

## Bangladesh shoudl be considered non-westernized
finalData$non_westernized[finalData$country == "BGD"] <- "yes"


finalData <- finalData %>%
  mutate(
    disease_group = case_when(
      disease == "healthy" ~ "Healthy",

      # Gut / digestion / liver / metabolic
      disease %in% c("acute_diarrhoea","adenoma","ascites;cirrhosis",
                     "ascites;cirrhosis;hepatitis","ascites;cirrhosis;hepatitis;schistosoma",
                     "ascites;cirrhosis;schistosoma","CDI","cirrhosis","cirrhosis;hepatitis",
                     "CMV;coeliac;gestational_diabetes","cystitis","fatty_liver",
                     "generic_diabetes","gestational_diabetes","gestational_diabetes;pre-eclampsia",
                     "hepatitis","IBD","IGT","infectiousgastroenteritis","pyelonefritis","pyelonephritis",
                     "salmonellosis", "STH","stomatitis","T1D","T1D;coeliac;irritable_bowel",
                     "T2D","tonsillitis","CRC","CRC;T2D") ~ "Gut / digestion",

      # Cardiovascular & Cancer
      disease %in% c("ACVD",
                     "melanoma","melanoma;metastases","melanoma;metastases;melanoma_surgery",
                     "melanoma;metastases_bone","melanoma;metastases_liver","melanoma;metastases_lung",
                     "melanoma;metastases_lung,metastases_nodes","melanoma;metastases_lung;metastases_adrenal",
                     "melanoma;metastases_lung;metastases_liver","melanoma;metastases_lung;metastases_liver;metastases_bone",
                     "melanoma;metastases_lung;metastases_liver;metastases_nodes",
                     "melanoma;metastases_lung;metastases_nodes",
                     "melanoma;metastases_lung;metastases_nodes;metastases_SQ",
                     "melanoma;metastases_nodes","melanoma;metastases_nodes;metastases_bone",
                     "melanoma;metastases_nodes;metastases_SQ","melanoma;metastases_SQ",
                     "melanoma;metastases_SQ;metastases_adrenal",
                     "melanoma;treatment_colitis",
                     "hypertension", "pre-eclampsia") ~ "Cardiovascular & Cancer",

      # Other / infections / immune / misc
      disease %in% c("asthma","BD","bronchitis","chorioamnionitis","cough","fever","ME/CFS",
                     "migraine","migraine;asthma","migraine;generic_diabetes","otitis","pneumonia",
                     "premature_born","RA","respiratoryinf","schizofrenia","sepsis","skininf","suspinf",
                     "MDRB") ~ "Other / infections",

      disease %in% "STEC" ~ "STEC",

      TRUE ~ "manual_check"
    )
  )

table(finalData$disease_group)

plot_heatmap <- function(df, anno_cols, num_prefix = c("hits_", "bt_"), ...) {
  num <- df %>%
    select(any_of(unlist(map(num_prefix,
                             ~ grep(.x, names(df), value = TRUE))))) %>%
    mutate(across(everything(), ~log10(.x + 1))) %>%
    as.data.frame()
  rownames(num) <- df$my_sample_id # assumes you have that column

  anno <- df %>%
    select(all_of(c("my_sample_id", anno_cols))) %>%
    as.data.frame()
  rownames(anno) <- anno$my_sample_id
  anno$my_sample_id <- NULL

  pheatmap::pheatmap(num,
                     annotation_row = anno,
                     show_rownames = FALSE,
                     ...)
}

png("figures/Loman2013_heat.png", res = 600, width = 8, height = 6,
    units = "in")
finalData %>%
  filter(DE2011 == "YES") %>%
  plot_heatmap(
    num_prefix = c("hits_.*stx2", "bt_.*stx2"),
    anno_cols = c("country", "disease", "stec_count", "shigatoxin_2_elisa",
                  "stool_texture" ),
    labels_col = c("Diamond 75% STX2A", "Diamond 75% STX2B",
                   "Diamond 95% STX2A", "Diamond 95% STX2B",
                   "Bowtie STX2A", "Bowtie STX2B")
    )
dev.off()

##
foo <- finalData %>%
    filter(DE2011 == "YES")

## In the newest iteration we are missing two positives!!! Why oh why? TODO!!
(nrow(foo[foo$STX2bt, ])+2)/(nrow(foo))

## Or did we process three more unifected samples? ... TODO!!!
nrow(foo[foo$STX2bt, ])/(nrow(foo)-3)

table(finalData[!finalData$study_name%in%"LomanNJ_2013", "STX2bt"], useNA="ifany")
### HERE DETECtion NEG = 16362    POS = 57

nrow(finalData[finalData$STX2bt & !finalData$study_name%in%"LomanNJ_2013", ])/
  nrow(finalData[!finalData$study_name%in%"LomanNJ_2013", ])*100
### 0.35 %

depth_model <- glm(STX2bt ~ number_bases,
                   family = binomial,
                   data = finalData[finalData$disease%in%"healthy", ])

summary(depth_model)


## From LomanNJ_2013 depth at their 0.67% sensitivity
# 1. Get reference prediction
ref_depth <- data.frame(number_bases = 1.075e9)
current_prob <- predict(depth_model, newdata = ref_depth, type = "response")

target_prob <- 0.67 ## TODO get it REAL
current_linear_pred <- predict(depth_model, newdata = ref_depth, type = "response")  # Log-odds
beta_1 <- coef(depth_model)["number_bases"]
new_intercept <- log(target_prob / (1 - target_prob)) - beta_1 * ref_depth$number_bases

depth_model_adjusted <- depth_model
depth_model_adjusted$coefficients["(Intercept)"] <- new_intercept

predict(depth_model_adjusted, newdata = ref_depth, type = "response")  # Should output 0.67

zero_depth <- data.frame(number_bases = 0)
predict(depth_model_adjusted, newdata = zero_depth, type = "response")  # Should output ZERO

new_data <- data.frame(number_bases = seq(1e8, 4e10, length.out = 100))
pred <- predict(depth_model_adjusted, newdata = new_data, type = "link", se.fit = TRUE)

new_data <- new_data %>%
  mutate(
    fit = pred$fit,                  # linear predictor (logit scale)
    se  = pred$se.fit,
    prob = plogis(fit),              # back-transform mean
    LL   = plogis(fit - 1.96 * se),  # CI lower bound
    UL   = plogis(fit + 1.96 * se)   # CI upper bound
  )

finalData[finalData$disease%in%"healthy", "pred_sensitivity"] <-
  predict(depth_model_adjusted,
          newdata = finalData[finalData$disease%in%"healthy","number_bases"],
          type="response")


# 4. Plot
extrapolation_plot <-
  ggplot(new_data, aes(x = number_bases, y = prob)) +
  geom_line() +
  geom_vline(xintercept = 1.075e9, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.67, linetype = "dashed", color = "red") +
  scale_x_continuous(
    limits = c(0, 4e10),
    labels = function(x) x / 1e9,
    breaks = seq(0, 4e10, by = 1e10)
  ) +
  ylim(0.6, 1) +
  theme_minimal() +
  labs(x = "Sequencing Depth (Gigabases)", y = "Adjusted Detection Probability",
       title = "Logistic model calibrated to 67% Sensitivity at 1.075Gb sequencing depth")

ggsave("figures/detection_depht.png", extrapolation_plot, bg = "white")

extrapolation_plot_studies <-
  extrapolation_plot +
  geom_jitter(
    data = finalData, alpha = 0.2, width = 1000000, height = 0.01,
    aes(x = number_bases, y = pred_sensitivity, color = study_name)
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1))  # force full opacity in legend
  ) +
  labs(
    color = "Study short name"   # custom legend title
  ) +
  theme(
    legend.position = c(0.98, 0.26),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6)),
    legend.key.size = unit(0.4, "cm"),   # smaller boxes
    legend.text = element_text(size = 7), # smaller labels
    legend.title = element_text(size = 8)
  )


ggsave("figures/detection_depht_studies.png", extrapolation_plot, bg = "white")


###  A feeling for how much data is in an area
# Center of ellipse
x0 <- 0.8e10   # pick around where prob ≈ 0.85 on your curve
y0 <- 0.85
a <- 6e9       # horizontal radius (adjust as needed)
b <- 0.15      # vertical radius

ellipse <- data.frame(
  x = x0 + a * cos(seq(0, 2*pi, length.out = 200)),
  y = y0 + b * sin(seq(0, 2*pi, length.out = 200))
)

inside <- finalData %>%
  filter(!is.na(number_bases), !is.na(pred_sensitivity)) %>%
  mutate(in_ellipse = ((number_bases - x0)^2 / a^2) +
           ((pred_sensitivity - y0)^2 / b^2) <= 1)

sum(inside$in_ellipse)   # count points

sum(inside$in_ellipse) / length(inside$in_ellipse)
## proportion of samples

sum(inside$number_bases[inside$in_ellipse]) ## count bases


sum(inside$number_bases[inside$in_ellipse])/
  sum(inside$number_bases) *100 ## proportion bases

extrapolation_plot_studies_ellipse <-
  extrapolation_plot_studies +
  geom_path(data = ellipse, aes(x = x, y = y), color = "black", linetype = "dashed")

ggsave("figures/detection_depht_studies_ellipse.png", extrapolation_plot_studies_ellipse,
       bg = "white")

extrapolation_plot_confidence <-
  extrapolation_plot +
  geom_ribbon(
    aes(x = number_bases, ymin = LL, ymax = UL),
    inherit.aes = FALSE,
    data = new_data,
    fill = "grey80", alpha = 0.5
  )

ggsave("figures/detection_depht_studies.png", extrapolation_plot_confidence,
       bg = "white")


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


disease_model<- glm(STX2bt ~ disease_group,
                    data=finalData,
                    family = binomial)

summary(disease_model)


west_model<- glm(STX2bt ~ non_westernized,
                    data=subset(finalData, !disease%in%"STEC"),
                    family = binomial)

summary(west_model)

table(subset(finalData, non_westernized == "yes")$disease,
      subset(finalData, non_westernized == "yes")$country)

tapply(finalData$age, finalData$non_westernized, fivenum, na.rm=TRUE)

country_model<- glm(STX2bt ~ country,
                 data=subset(finalData, !disease%in%"STEC"),
                 family = binomial)

summary(country_model)

finalData$stunting <- finalData$disease%in%"STH"

stunting_model <- glm(STX2bt ~ stunting,
                       data=subset(finalData, !disease%in%"STEC" &
                                     finalData$non_westernized == "yes"
                       ),
                       family = binomial)

### nothing for stunting not even within non-westernized

finalData %>%
  dplyr::filter(DE2011 == "NO" &
                STX2bt == TRUE) %>%
  plot_heatmap(anno_cols = c("non_westernized", "disease_group"))

# Summarise counts AND non_westernized together
country_summary <- finalData %>%
  group_by(country) %>%
  summarise(
    n_samples = n(),
    non_westernized = any(non_westernized == "yes"),
    .groups = "drop"
  )

# Join to world map
world_data <- world %>%
  left_join(country_summary, by = c("iso_a3" = "country"))

# 4. Plot with ggplot2
country_plot <-
  ggplot(world_data) +
  geom_sf(aes(fill = n_samples), color = "grey80") +
  scale_fill_viridis_c(option = "plasma", trans = "log", na.value = "grey90",
                       breaks = c(1, 10, 100, 1000, 10000)
  ) +
  ylim(-60, 80) +
  theme_minimal() +
  labs(fill = "Samples",
       title = "Number of samples per country",
       subtitle = "log scale color gradient")

ggsave("figures/country_plot.png", country_plot, bg = "white", width = 9,
       height = 4, dpi = 600)

country_plot_western <-
  country_plot +
  geom_sf(
  data = world_data %>% filter(non_westernized == TRUE),
  fill = NA, color = "black", size = 1) +
  labs(fill = "Samples",
       title = "Number of samples per country",
       subtitle = "highglihgting \"non westernized countries\"")

ggsave("figures/country_plot_western.png", country_plot_western,
       bg = "white", width = 9, height = 4, dpi = 600)
