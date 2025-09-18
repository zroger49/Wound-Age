#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Parse the metadata and perform exploration analysis
# @software version: R=4.4.2

# Outputs the final version of the metadata

# Load libraries ----
spsm <- suppressPackageStartupMessages
spsm(library(tidyverse))
spsm(library(data.table))
spsm(library(reshape2))
spsm(library(ggcorrplot))
spsm(library(ComplexHeatmap))
spsm(library(circlize))
spsm(library(RColorBrewer))
spsm(library(caret))
spsm(library(lme4))

# Create output dir ----
out_data <- "results_new/01_data_exploration"
if (!dir.exists(out_data)){
  dir.create(out_data, recursive  =T )
}

# Load data ----
# (RPM). RPM is the reads normalized the total reads in the sample
rpm1 <- fread("data/rpm/Auto_IvaGomes_THGE_29-05-2024_B_Chip1_torrent-server_6371.rpm.bcmatrix.txt") %>% column_to_rownames("Gene")
rpm2 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_A_Chip2_torrent-server_21.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")
rpm3 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_B_Chip3_torrent-server_23.rpm.bcmatrix.xls.txt")%>% column_to_rownames("Gene")
rpm4 <- fread("data/rpm/Auto_IvaGomes_THGE_01-07-2024_D_Chip4_torrent-server_46.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")
#Modify rpm4 names
colnames(rpm4)[2:5] <- paste0(colnames(rpm4)[2:5], "2")

a <- list(rpm1, rpm2, rpm3, rpm4)
rpms <- lapply(a, function(x) x %>% select(-c(Target, COSMIC_CGC_FLAG, NCBI_NAME, HGNC_SYMBOL_ACC, MIM_MORBID_DESCRIPTION, ENTREZ_GENE_ID, U133PLUS2_PSID)))
rpms <- do.call("cbind", rpms)

# Remove two duplicated samples 

rpms <- rpms %>% select(-c( IonXpress_037, IonXpress_038))

# Load matching samples to replace colnames 
match <- list(fread("data/expression_data/chip1_matching_samples.txt"), 
              fread("data/expression_data/chip2_matching_samples.txt"), 
              fread("data/expression_data/chip3_matching_samples.txt"), 
              fread("data/expression_data/chip4_matching_samples_mod.txt")
)

match <- read.csv("data/match.csv") %>% select(-Sample.Name)
match_vector <- match$Sample_name
names(match_vector) <- match$Barcode.ID

# Replace the names in the expression table
colnames(rpms) <- match_vector[colnames(rpms)]
# Load metadata
metadata_sample.full <- read.csv("data/metadata.csv", row.names = 1, sep = ";") 
metadata_sample.full <- metadata_sample.full[-c(21,22),]


# Basic transforms and selection ----

# location mapping
location_mapping <- function(raw_locations) {
  normalized <- tolower(trimws(raw_locations))
  dplyr::case_when(
    grepl("chest|thoraco|upper chest", normalized) ~ "Chest",
    grepl("abdomen|body center|flank|groin|hip|pelvis", normalized) ~ "Abdomen/Pelvis",
    grepl("arm|hand", normalized) ~ "Upper Limb",
    grepl("leg|thigh|tight|femoral", normalized) ~ "Lower Limb",
    grepl("neck", normalized) ~ "Neck",
    grepl("\\?|or", normalized) ~ "Ambiguous",
    TRUE ~ "Unknown"
  )
}

metadata_sample.full$body_location <- location_mapping(metadata_sample.full$localisation)
metadata_sample.full$Chip <- as.character(metadata_sample.full$Chip)

metadata_sample.full.selected <- metadata_sample.full %>% 
  mutate(group = ifelse(substr(metadata_sample.full$Sample, 1, 1) == "R", "Control", "Wound")) %>%
  select(c(Subject, Sample_name, Age, Gender, RIN, Chip, PMI, body_location, case.type,group))

categorical_vars <- metadata_sample.full.selected %>%
  select(where(~ is.character(.x) || is.factor(.x))) %>%
  names()
categorical_vars <- categorical_vars[!categorical_vars %in% c("Subject", "Sample_name")]

# Replace empty strings by NA
metadata_sample.full.selected <- metadata_sample.full.selected %>% as.data.frame() %>%
  mutate(across(where(~ is.character(.) | is.factor(.)), ~ na_if(.x, "")))

encoded_df <- metadata_sample.full.selected

# One-hot encode categorical variables with NAs preserved
for (var in categorical_vars) {
  dummies <- dummyVars(reformulate(var), data = metadata_sample.full.selected)
  one_hot_data <- predict(dummies, newdata = metadata_sample.full.selected)
  colnames(one_hot_data) <- gsub("_", " ", colnames(one_hot_data))
  
  # Mark rows that were all-zero => original NA; replace with NA rows
  NA_rows <- rowSums(one_hot_data) == 0
  if (any(NA_rows)) one_hot_data[NA_rows, 1:ncol(one_hot_data)] <- NA
  
  if (ncol(one_hot_data) == 2){
    one_hot_data <- one_hot_data[,1, drop = FALSE]
  }
  
  encoded_df <- encoded_df %>%
    select(-all_of(var)) %>%
    bind_cols(as.data.frame(one_hot_data))
}

# Unsupervised analysis ----
## Hierarchical clustering (Wound samples only)

# Log transform full rpm matrix (used later for PCA on all samples)
rpms_log <- log(rpms + 1)

# Subset to wound samples only
encoded_df.wound <- encoded_df[encoded_df$groupControl == "0", ]
encoded_df.wound <- encoded_df.wound[encoded_df.wound$Sample %in% colnames(rpms_log), ]

# Prepare rpms for wound samples (log version)
rpms.log.wound <- rpms_log %>% select(all_of(encoded_df.wound$Sample))

# Create a wound age factor, dividing the current samples
metadata_sample.wound_age <- metadata_sample.full %>% select(Subject,  wound.age) %>% distinct()
metadata_sample.wound_age <- metadata_sample.wound_age %>% rename(wound_age = `wound.age`)

metadata_sample.wound_age$wound_age <- c("2/3 days", "2/3 days",
                                         "7 days", "12-13 days", 
                                         "minutes-hours","minutes-hours", 
                                         "2/3 days", "minutes-hours", 
                                         "1 day", ">20 days", 
                                         "12-13 days", "2/3 days", 
                                         "minutes-hours", "minutes-hours", 
                                         "minutes-hours", "minutes-hours", 
                                         ">20 days","minutes-hours", 
                                         "minutes-hours", "minutes-hours", 
                                         "minutes-hours","2/3 days", 
                                         "7 days")

metadata_sample.wound_age$wound_age <- factor(
  metadata_sample.wound_age$wound_age,
  levels = c("minutes-hours", "1 day", "2/3 days", "7 days", "12-13 days", ">20 days"),
  ordered = TRUE
)

# Merge wound_age into sample-level metadata for wound samples
metadata_sample.full.selected.wound <- metadata_sample.full.selected %>%
  filter(group == "Wound") %>%
  left_join(metadata_sample.wound_age, by = "Subject")


## Hierarchical clustering  ----
## Remove the effect of other variables with a linear regression

# Compute residuals per gene using lmer (preserve random effect Chip)
residuals_mat <- apply(rpms.log.wound, 1, function(x){
  all_data <- metadata_sample.full.selected.wound
  all_data$expression <- x
  resid(lmer(expression ~ Age + Gender + RIN + PMI + case.type + (1|Chip), data = all_data))
})
residuals_mat <- t(residuals_mat)
saveRDS(residuals_mat, file = file.path(out_data, "residuals_model.rds"))

# select top 5000 variable genes from residuals
sd_genes_res <- apply(residuals_mat, 1, sd)
top5000_var_genes_res <- order(sd_genes_res, decreasing = TRUE)[1:5000]
residuals.top <- residuals_mat[top5000_var_genes_res, ] %>% t()
residuals.top.scaled <- scale(residuals.top, center = TRUE, scale = TRUE)
row.names(residuals.top) <- colnames(rpms.log.wound)

dist_mat_res <- dist(residuals.top.scaled)
hc_res <- hclust(dist_mat_res, method = "ward.D2")


# Colors / annotation for heatmap 
gender_colour <-  brewer.pal(n = 3, "Set1")[c(1,3)]
chip_colour <- brewer.pal(n = 10, "Set3")[c(1,3,5,7)]
body_section_cols <- c("#E41A1C","#377EB8","#4DAF4A", "#984EA3","#FF7F00")
wound_time_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

names(gender_colour) <- unique(metadata_sample.full.selected.wound$Gender)
names(chip_colour) <- unique(metadata_sample.full.selected.wound$Chip)
names(body_section_cols) <- unique(metadata_sample.full.selected.wound$body_location)
names(wound_time_cols) <- unique(metadata_sample.full.selected.wound$wound_age)

ha_res <- HeatmapAnnotation(
  `Wound time` = metadata_sample.full.selected.wound$wound_age,
  `Body section` = metadata_sample.full.selected.wound$body_location,
  Sex = metadata_sample.full.selected.wound$Gender,
  Chip = metadata_sample.full.selected.wound$Chip,
  col = list(
    `Body section` = body_section_cols,
    Sex = gender_colour,
    Chip = chip_colour,
    `Wound time` = wound_time_cols
  ),
  annotation_name_side = "left"
)

row.names(residuals.top.scaled) <- substr(row.names(residuals.top), 2, 4)

ht_res <- Heatmap(
  t(residuals.top.scaled),
  cluster_columns = as.dendrogram(hc_res),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = T,
  top_annotation = ha_res,
  column_title = "HC on residuals (top 5000 variable genes)",
  heatmap_legend_param = list(title = "Residual Expression"),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

pdf(file.path(out_data, "01_HC_residuals_top5000_genes.pdf"), w = 12)
draw(ht_res)
dev.off()

# Save residuals-based HC data
saveRDS(list(
  residuals = residuals_mat,
  residuals_top = residuals.top,
  residuals_top_scaled = residuals.top.scaled,
  hc_res = hc_res
), file = file.path(out_data, "hc_residuals_inputs_outputs.rds"))

# Assign a group based on the clusters
metadata_sample.wound_age <- metadata_sample.wound_age %>%
  mutate(
    wound_time = case_when(
      wound_age == "minutes-hours"             ~ "Acute",
      wound_age %in% c("1 day", "2/3 days")   ~ "Intermediary",
      wound_age %in% c("7 days", "12-13 days", ">20 days") ~ "Prolonged",
      TRUE ~ NA_character_
    )
  )

# Add time group to the enconded variable 
dummies <- dummyVars(~ wound_time, data = metadata_sample.wound_age)
wound_time_encoded <- predict(dummies, newdata = metadata_sample.wound_age) %>% as.data.frame()
colnames(wound_time_encoded) <- gsub("_", " ", colnames(wound_time_encoded))

wound_time_encoded$Subject <- metadata_sample.wound_age$Subject


## PCA----
## Keep top genes across all samples
sd_genes_all <- apply(rpms_log, 1, sd)
top5000_var_genes_all <- order(sd_genes_all, decreasing = TRUE)[1:5000]
rpms_top <- rpms_log[top5000_var_genes_all, ] %>% t()
rpms_top <- scale(rpms_top, center = TRUE, scale = TRUE)

# Run PCA
pca_data <- prcomp(rpms_top)
pc_scores <- pca_data$x[, 1:10]
pc_var <- pca_data$sdev^2
pc_var_explained <- round(100 * pc_var / sum(pc_var), 1)  # % variance
pc_names <- paste0("PC", 1:10, " (", pc_var_explained[1:10], "%)")
colnames(pc_scores) <- pc_names
rownames(pc_scores) <- rownames(rpms_top)

# Add wound time information
metadata_aligned <- encoded_df %>% 
  merge(wound_time_encoded, by = "Subject")

metadata_aligned <-  metadata_aligned[match(rownames(pc_scores), metadata_aligned$Sample_name), ]

# Align metadata to PCA scores
pc_meta <- cbind(pc_scores, metadata_aligned)

# remove unnecessary columns
metadata_aligned <- metadata_aligned %>% 
  select(-c(Subject, Sample_name, `body locationChest`, `body locationLower Limb`, `body locationNeck`, `body locationUpper Limb`, `body locationAbdomen/Pelvis`))

# Rename metadata columns (hardcoded as before)
colnames(metadata_aligned) <- c("Age", "RIN", "PMI", "Gender Female", "CHIP 1", "CHIP 2", "CHIP 3", "CHIP 4", "case criminal", "Group control", "Wound acute", "Wound intermediate", "Wound prolonged")
metadata_aligned <- metadata_aligned[,c(1:9,11:13,10)]



# compute correlations and p-values
cor_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_matrix) <- colnames(pc_scores)
colnames(cor_matrix) <- colnames(metadata_aligned)

cor_pvalues_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_pvalues_matrix) <- colnames(pc_scores)
colnames(cor_pvalues_matrix) <- colnames(metadata_aligned)

for (i in 1:ncol(pc_scores)) {
  for (j in 1:ncol(metadata_aligned)) {
    x <- pc_scores[,i]
    y <- metadata_aligned[,j]
    cor_matrix[i, j] <- cor(x, y, use = "complete.obs", method = "spearman")
    cor_pvalues_matrix[i,j] <- cor.test(x, y, use = "complete.obs", method = "spearman")$p.value
  }
}

# adjust p-values across all tests (FDR)
cor_pvalues_adj_matrix <- matrix(
  p.adjust(as.vector(cor_pvalues_matrix), method = "fdr"),
  nrow = nrow(cor_pvalues_matrix),
  ncol = ncol(cor_pvalues_matrix)
)
rownames(cor_pvalues_adj_matrix) <- rownames(cor_pvalues_matrix)
colnames(cor_pvalues_adj_matrix) <- colnames(cor_pvalues_matrix)

# create data.frame for plotting
cor_df <- melt(cor_matrix, varnames = c("PC", "Metadata"), value.name = "cor")
pvalue_df <- melt(cor_pvalues_adj_matrix, varnames = c("PC", "Metadata"), value.name = "FDR")
cor_df$FDR <- pvalue_df$FDR

cor_df$Metadata <- factor(cor_df$Metadata, levels = c("Age", "Gender Female", "RIN", "PMI", "case criminal", 
                                                      "CHIP 1", "CHIP 2", "CHIP 3", "CHIP 4",
                                                      "Wound control", "Wound acute", 
                                                      "Wound intermediate", "Wound prolonged",
                                                      "Group control"))

pdf(file.path(out_data, "02_PC_vs_metadata.pdf"), w = 10, h = 8)
ggplot(cor_df, aes(x = Metadata, y = PC, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Correlation / R", limits = c(-1, 1)) +
  theme_minimal() +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Association between First 10 PCs and Metadata") + 
  xlab("") + ylab("")
dev.off()

pdf(file.path(out_data, "02_PC_vs_metadat_clean.pdf"), w = 8, h = 6)
ggplot(cor_df, aes(x = Metadata, y = PC, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(FDR < 0.05, sprintf("%.2f", cor), "")), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Correlation", limits = c(-1, 1)) +
  theme_minimal() +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Association between First 10 PCs and Metadata (FDR < 0.05 shown)") + 
  xlab("") + ylab("")  + 
  theme(axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"))
dev.off()

# Save data
saveRDS(list(
  pca_data = pca_data,
  pc_scores = pc_scores,
  pc_var_explained = pc_var_explained,
  cor_df = cor_df,
  cor_matrix = cor_matrix,
  cor_pvalues_adj_matrix = cor_pvalues_adj_matrix,
  metadata_aligned = metadata_aligned
), file = file.path(out_data, "pca_and_correlations.rds"))

## Plot the PCA ----
# Make sure metadata aligns with PCA scores
metadata_for_plot <- metadata_sample.full[match(rownames(pc_scores), metadata_sample.full$Sample), ]

metadata_for_plot$wound_time <- factor(metadata_for_plot$wound_time, 
                                       levels = c("control","Acute", "Intermediary", "Prolonged"))


pca_plot_scores <- pc_scores
colnames(pca_plot_scores)[1:2] <- c("PC1", "PC2")

# PCA plot: PC1 vs PC2
pca_plot <- ggplot(data = pca_plot_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = metadata_for_plot$wound_time), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("control" = "#F0E442",
                                "Acute" = "#E69F00",
                                "Intermediary" = "#56B4E9",
                                "Prolonged" = "#009E73")) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of Top 5000 Variable Genes",
       x = paste0("PC1 (", pc_var_explained[1], "% variance)"),
       y = paste0("PC2 (", pc_var_explained[2], "% variance)"),
       color = "Wound Time") +
  theme(legend.position = "right")

pdf(file.path(out_data, "02_PCA.pdf"), w = 8, h = 6)
print(pca_plot)
dev.off()

# Metadata exploration ----

#Sample distribution by clinical / police case

metadata_sample.full <- metadata_sample.full %>% 
  merge(metadata_sample.wound_age, by.x ="Subject", by.y = "Subject") %>% 
  mutate(wound_time = ifelse(substr(Sample_name,1,1) == "W", wound_time, "control"))

write.csv(metadata_sample.full, "data/metadata_wound_time_assigned.csv")
#wound_time_encoded

plot <- metadata_sample.full %>% 
  ggplot(aes(x = wound_time)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribution of Wound group",
    x = "",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~case.type)


pdf(paste0(out_data, "/03.wound_time_per_casetype.pdf"), w = 6, h = 3)
print(plot)
dev.off()

## Make a plot of the correlation between variables
encoded_df_to_cor <- encoded_df %>% select(-c(Subject, Sample_name))

cor_mat <- cor(metadata_aligned, use = "pairwise.complete.obs", method = "spearman")


#row.names(df_encoded) <- row.names(metadata)
pdf(file.path(out_data, "04_corr_metadata.pdf"), w = 9, h = 10)
ggcorrplot(cor_mat, 
           method = "square", 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           colors = c("blue", "white", "red"),
           title = "Correlation Heatmap (One-Hot Encoded)",
           ggtheme = theme_minimal())

dev.off()


# Plot sample size ----

# Make sure wound_time is a factor with levels in order
metadata_sample.full$wound_time <- factor(metadata_sample.full$wound_time,
                                          levels = c("control", "Acute", "Intermediary", "Prolonged"))

# Create a grouping variable: Control vs Wounds
metadata_sample.full <- metadata_sample.full %>%
  mutate(group_bar = ifelse(wound_time == "control", "Control", "Wound"))

# Count samples per wound_time per group
plot_data <- metadata_sample.full %>%
  group_by(group_bar, wound_time) %>%
  summarise(count = n(), .groups = "drop")

grey_palette <- c(
  "control"     = "#D9D9D9", 
  "Acute"       = "#9E9E9E",  
  "Intermediary"= "#6E6E6E",  
  "Prolonged"   = "#3B3B3B"   
)


# Make the barplot
pdf(file.path(out_data, "sample_size.pdf"), w = 4, h = 4)
ggplot(plot_data, aes(x = group_bar, y = count, fill = wound_time)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5),  # center of each stacked bar
            size = 5, colour = "black") +
  scale_fill_manual(values = grey_palette) +
  labs(x = "", y = "Number of samples", fill = "Wound Time") +
  theme_classic()+
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        legend.position = "right")
dev.off()
