#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Find Biomarkers
# @software version: R=4.4.2

# Library
sh <- suppressPackageStartupMessages
sh(library(tidyverse))
library(data.table)

#Set res 
out_data <- "results_new/06_biomarkers"
if (!dir.exists(out_data)){
  dir.create(out_data)
}


# Load the data
# Load data (RPM). RPM is the reads normalized the total reads in the sample 
rpm1 <- fread("data/rpm/Auto_IvaGomes_THGE_29-05-2024_B_Chip1_torrent-server_6371.rpm.bcmatrix.txt") %>% column_to_rownames("Gene")
rpm2 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_A_Chip2_torrent-server_21.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")
rpm3 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_B_Chip3_torrent-server_23.rpm.bcmatrix.xls.txt")%>% column_to_rownames("Gene")
rpm4 <- fread("data/rpm/Auto_IvaGomes_THGE_01-07-2024_D_Chip4_torrent-server_46.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")

#Modify rpm4 names
colnames(rpm4)[2:5] <- paste0(colnames(rpm4)[2:5], "2")

a <- list(rpm1, rpm2, rpm3, rpm4)
rpms <- lapply(a, function(x) x %>% select(-c(Target, COSMIC_CGC_FLAG, NCBI_NAME, HGNC_SYMBOL_ACC, MIM_MORBID_DESCRIPTION, ENTREZ_GENE_ID, U133PLUS2_PSID)))
rpms <- do.call("cbind", rpms)

# Load matching samples to replace colnames 
match <- list(fread("data/expression_data/chip1_matching_samples.txt"), 
              fread("data/expression_data/chip2_matching_samples.txt"), 
              fread("data/expression_data/chip3_matching_samples.txt"), 
              fread("data/expression_data/chip4_matching_samples_mod.txt")
)

match <- do.call("rbind", match)
match_vector <- match$`Sample Name`
names(match_vector) <- match$`Barcode ID`

# Replace the names in the expression table
colnames(rpms) <- match_vector[colnames(rpms)]
#any(duplicated(colnames(rpms)))


# Load metadata 
metadata_subject <- fread("data/metadata_full.csv") %>% select(-c(V23, V24))
metadata_sample <- fread("data/metadata_complete.csv")

subjects_with_seq_data <-  substring(colnames(rpms), 2)
subjects_with_seq_data <- gsub("\\.\\d", "", subjects_with_seq_data) # Removes ".1" 
metadata_subject.seq_samples <- metadata_subject %>% filter(`Wound sample` %in% subjects_with_seq_data) 
dim(metadata_subject.seq_samples) # [1] 23 22
# Merge with the sample level metadata 

metadata_sample$Subject <- substring(metadata_sample$Sample, 2)
metadata_sample$Subject <- gsub("\\.\\d", "", metadata_sample$Subject) # Removes ".1" 

metadata_sample.full <- merge(metadata_sample, metadata_subject.seq_samples, by.x = "Subject", by.y = "Wound sample")

# Manually add the early, intermediary and late
metadata_sample.full$timepoint_analysis <- rep(c("Intermediary","Intermediary", "Late", "Late", "Early", "Early", 
                                                 "Intermediary", "Early", "Intermediary", "Late", "Late","Late","Intermediary", 
                                                 "Early", "Early", "Early", "Early", "Late", "Early", "Early", "Early", "Early", 
                                                 "Intermediary", "Late"),each = 2)

#Load DEGs 
degs_file <- c(
  "results_new/02_degs/control_late.txt",
  "results_new/02_degs/control_early.txt",
  "results_new/02_degs/control_intermediary.txt"
)

DEGS <- lapply(degs_file, function(x) read.csv(x, sep = "\t", skip = 4, dec = ","))
DEGs_control_late <- DEGS[[1]]
DEGs_control_early <- DEGS[[2]]
DEGs_control_intermediary <- DEGS[[3]]

# Ratio of expression
rpms.w <- rpms[,c(1:ncol(rpms)) %% 2 == 0]
rpms.c <- rpms[,c(1:ncol(rpms)) %% 2 == 1]

ratio.log <- log2((rpms.w + 0.1) / (rpms.c +0.1))
ratio <- (rpms.w + 0.1) / (rpms.c +0.1)


#Differential expressed exclusively in that group
#Ratio >= 2 in all the sample (of that group). 
#Basal expression in the wound high or Basal expression in Wound high (positive and negative markers) 

#In all, I tried to get ~20 genes per group. 

degs_early <- DEGs_control_early %>% filter(FDR.P.val < 0.05)
degs_int <- DEGs_control_intermediary %>% filter(FDR.P.val < 0.05)
degs_late <- DEGs_control_late %>% filter(FDR.P.val < 0.05)


### Early
#### Get exclusive genes
degs_early_exclusive <- setdiff(degs_early$ID, union(degs_int$ID, degs_late$ID))
degs_early_exclusive <- degs_early %>% filter(ID %in% degs_early_exclusive)

### Get ratio for those genes 
#### I am mostly interested in positive markers, so no ploting in Log is fine here.
metadata_sample.full.early <- metadata_sample.full %>% filter(timepoint_analysis == "Early")
early_samples <- metadata_sample.full.early %>% pull(Sample)
early_samples <- early_samples[grepl("^W", early_samples)]

ratio_early_genes <- ratio[degs_early_exclusive$ID,early_samples]

###Median ratios 
median_early_degs_ratio <-  apply(ratio_early_genes, 1, median)

plot(hist(median_early_degs_ratio, 100)) # Lot of genes with a high effect, without taking into account confounding factors. 

# Filter for genes with a high effect
degs_early_exclusive.high.effect <- degs_early_exclusive %>% filter(abs(Fold.Change) > 2)
rpms.w.early.genes.higheffect <- rpms.w[degs_early_exclusive.high.effect$ID, early_samples]
# Filter some lowexpressed genes
apply(rpms.w.early.genes.higheffect,1,median)[order(apply(rpms.w.early.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.early.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])

degs_early_exclusive.high.expression <- degs_early_exclusive.high.effect %>% filter(ID %in% keep)

# Ratio for these genes
ratio.early.candidate.markers <- ratio_early_genes[degs_early_exclusive.high.expression$ID,]
# order by median ratio
ratio.early.candidate.markers <- ratio.early.candidate.markers[order(apply(ratio.early.candidate.markers, 1, median), decreasing = T), ]

# Number of samples with ratio >2 
good_samples <- rowSums(ratio.early.candidate.markers > 2)

data_sorted <- sort(good_samples, decreasing = TRUE)

pdf(file.path(out_data, "01_early_biomarkers_sample_number.pdf"), w = 8, h = 6)
barplot(
  data_sorted,
  col = "skyblue",
  border = "black",
  las = 2, # rotate labels
  xlab = "",
  ylab = "Number of Samples with ratio > 2"
)
dev.off()

plot(good_samples, 
        degs_early_exclusive.high.expression$Fold.Change)

# Biomarkers to plot
ratio.early.markers <- ratio.early.candidate.markers[good_samples >= 8,] # Select genes with ratio in more that 8 samples


# Reshape to long format
ratio.early.candidate.markers.long <- ratio.early.markers %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "ratio"
  )

# Merge with metadata
ratio.early.candidate.markers.long <- ratio.early.candidate.markers.long %>%
  left_join(metadata_sample.full.early, by = c("sample" = "Sample"))


wound_age_order <- c(
  "seconds",
  "seconds/ minutes",
  "minutes", 
  "min", 
  "15 min",
  "minutes - < 1 hour",
  "hours (< 10h)", 
  "hours"
)

ratio.early.candidate.markers.long <- ratio.early.candidate.markers.long %>% 
  mutate(`wound age` = factor(`wound age`, levels = wound_age_order)) %>% 
  arrange(`wound age`) %>% 
  mutate(sample = factor(sample, levels = unique(sample)))
                                                    

plots <- lapply(unique(ratio.early.candidate.markers.long$gene), function(x){
  ratio.early.candidate.markers.long %>% 
    filter(gene == x) %>%
    ggplot(aes(x = sample, y = ratio)) +
    geom_point(size = 5, color = "steelblue") +
    geom_text(
      aes(label = paste(localisation, "\n", `Wound age`, age, Gender, sep = ", ")),
      vjust = -0.5, size = 2.5
    ) +
    labs(
      title = paste("Expression for", x),
      x = "",
      y = "Expression Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    #expand_limits(y = max(ratio) * 1.1) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 14),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 20) 
    ) + 
    coord_cartesian(clip = "off")
})

for (i in seq_along(plots)) {
  ggsave(file.path(out_data, paste0("02_plot_", unique(ratio.early.candidate.markers.long$gene)[i], ".early.pdf")), plots[[i]], width = 10, height = 6)
}



### Intermediary
#### Get exclusive genes
degs_int_exclusive <- setdiff(degs_int$ID, union(degs_early$ID, degs_late$ID))
degs_int_exclusive <- degs_int %>% filter(ID %in% degs_int_exclusive)

### Get ratio for those genes 
#### I am mostly interested in positive markers, so no ploting in Log is fine here.
metadata_sample.full.int <- metadata_sample.full %>% filter(timepoint_analysis == "Intermediary")
int_samples <- metadata_sample.full.int %>% pull(Sample)
int_samples <- int_samples[grepl("^W", int_samples)]

ratio_int_genes <- ratio[degs_int_exclusive$ID,int_samples]

###Median ratios 
median_int_degs_ratio <-  apply(ratio_int_genes, 1, median)

plot(hist(median_int_degs_ratio, 100)) # Lot of genes with a high effect, without taking into account confounding factors. 

# Filter for genes with a high effect
degs_int_exclusive.high.effect <- degs_int_exclusive %>% filter(abs(Fold.Change) > 2)
# This time I got both positive and negative markers.
degs_int_exclusive.high.effect.w <- degs_int_exclusive.high.effect %>% filter(Fold.Change < 0)
rpms.w.int.genes.higheffect <- rpms.w[degs_int_exclusive.high.effect.w$ID, int_samples]

degs_int_exclusive.high.effect.c <- degs_int_exclusive.high.effect %>% filter(Fold.Change > 0)
rpms.c.int.genes.higheffect <- rpms.c[degs_int_exclusive.high.effect.c$ID, ]


# Filter some lowexpressed genes
apply(rpms.w.int.genes.higheffect,1,median)[order(apply(rpms.w.int.genes.higheffect,1,median))]
apply(rpms.c.int.genes.higheffect,1,median)[order(apply(rpms.w.int.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.int.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])
degs_int_exclusive.high.expression <- degs_int_exclusive.high.effect %>% filter(ID %in% keep)

median_rpm <- apply(rpms.c.int.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])
degs_int_exclusive.high.expression <- rbind(
  degs_int_exclusive.high.expression,
  degs_int_exclusive.high.effect %>% filter(ID %in% keep)
)

degs_int_exclusive.high.expression.candidate.markers <- degs_int_exclusive.high.expression %>% 
  arrange(-abs(Fold.Change)) %>%
  slice_head(n = 30)

# Ratio for these genes
ratio.int.candidate.markers <- ratio_int_genes[degs_int_exclusive.high.expression.candidate.markers$ID,]
# order by median ratio
ratio.int.candidate.markers <- ratio.int.candidate.markers[order(apply(ratio.int.candidate.markers, 1, median), decreasing = T), ]

# Number of samples with ratio >2 
good_samples.pos <- rowSums(ratio.int.candidate.markers > 2)
good_sample.neg <- rowSums(ratio.int.candidate.markers < 0.5)


good_samples.pos <- good_samples.pos[good_samples.pos >= 5]
good_sample.neg <- good_sample.neg[good_sample.neg >= 5]

good_samples <- c(good_sample.neg, good_samples.pos)


# Biomarkers to plot
ratio.int.markers <- ratio.int.candidate.markers[names(good_samples),] # Select genes with ratio in more that 5 samples

# Reshape to long format
ratio.int.candidate.markers.long <- ratio.int.markers %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "ratio"
  )

# Merge with metadata
ratio.int.candidate.markers.long <- ratio.int.candidate.markers.long %>%
  left_join(metadata_sample.full.int, by = c("sample" = "Sample"))


wound_age_order <- c(
  "24 hours",
  "2 days",
  "55 hours (2,3 days)", 
  "70 hours (2,9 days)"
)

ratio.int.candidate.markers.long <- ratio.int.candidate.markers.long %>% 
  mutate(`wound age` = factor(`wound age`, levels = wound_age_order)) %>% 
  arrange(`wound age`) %>% 
  mutate(sample = factor(sample, levels = unique(sample)))


plots <- lapply(unique(ratio.int.candidate.markers.long$gene), function(x){
  ratio.int.candidate.markers.long %>% 
    filter(gene == x) %>%
    ggplot(aes(x = sample, y = ratio)) +
    geom_point(size = 5, color = "steelblue") +
    geom_text(
      aes(label = paste(localisation, "\n", `Wound age`, age, Gender, sep = ", ")),
      vjust = -0.5, size = 2.5
    ) +
    labs(
      title = paste("Expression for", x),
      x = "",
      y = "Expression Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    #expand_limits(y = max(ratio) * 1.1) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 14),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 20) 
    ) + 
    coord_cartesian(clip = "off")
})

for (i in seq_along(plots)) {
  ggsave(file.path(out_data, paste0("03_plot_", unique(ratio.int.candidate.markers.long$gene)[i], ".intermediary.pdf")), plots[[i]], width = 6, height = 6)
}

### Late
#### Get exclusive genes
degs_late_exclusive <- setdiff(degs_late$ID, union(degs_int$ID, degs_early$ID))
degs_late_exclusive <- degs_late %>% filter(ID %in% degs_late_exclusive)

### Get ratio for those genes 
#### I am mostly interested in positive markers, so no ploting in Log is fine here.
metadata_sample.full.late <- metadata_sample.full %>% filter(timepoint_analysis == "Late")
late_samples <- metadata_sample.full.late %>% pull(Sample)
late_samples <- late_samples[grepl("^W", late_samples)]

ratio_late_genes <- ratio[degs_late_exclusive$ID,late_samples]

###Median ratios 
median_late_degs_ratio <-  apply(ratio_late_genes, 1, median)

plot(hist(median_late_degs_ratio, 100)) # Lot of genes with a high effect, without taking into account confounding factors. 

# Filter for genes with a high effect
degs_late_exclusive.high.effect <- degs_late_exclusive %>% filter(abs(Fold.Change) > 2)
rpms.w.late.genes.higheffect <- rpms.w[degs_late_exclusive.high.effect$ID, late_samples]
#Most of these are positive markers, only 1 is a negative marker. Since the expression of that genes seems fine, I do not need to filter for it in separate
# Filter some lowexpressed genes
apply(rpms.w.late.genes.higheffect,1,median)[order(apply(rpms.w.late.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.late.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])

degs_late_exclusive.high.expression <- degs_late_exclusive.high.effect %>% filter(ID %in% keep)

# Ratio for these genes
ratio.late.candidate.markers <- ratio_late_genes[degs_late_exclusive.high.expression$ID,]
# order by median ratio
ratio.late.candidate.markers <- ratio.late.candidate.markers[order(apply(ratio.late.candidate.markers, 1, median), decreasing = T), ]


# Number of samples with ratio >2 
good_samples <- rowSums(ratio.late.candidate.markers > 2)
#AMOTL2 is the exception, because its a negative marker
data_sorted <- sort(good_samples, decreasing = TRUE)
data_sorted[14] <- 3 # Visual inspection of the gene reveals 3 samples where the expression is higher in control with a ratio of 0.5 or lower. 

pdf(file.path(out_data, "01_late_biomarkers_sample_number.pdf"), w = 8, h = 6)
barplot(
  data_sorted,
  col = "skyblue",
  border = "black",
  las = 2, # rotate labels
  xlab = "",
  ylab = "Number of Samples with ratio > 2"
)
dev.off()

plot(good_samples, 
     degs_late_exclusive.high.expression$Fold.Change)

# Biomarkers to plot
ratio.late.markers <- ratio.late.candidate.markers[good_samples >= 6,] # Select genes with ratio in more that 6 samples

# Reshape to long format
ratio.late.candidate.markers.long <- ratio.late.markers %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "ratio"
  )

# Merge with metadata
ratio.late.candidate.markers.long <- ratio.late.candidate.markers.long %>%
  left_join(metadata_sample.full.late, by = c("sample" = "Sample"))


wound_age_order <- c(
  "7 days",
  "12 days",
  "12/13 days", 
  "22 days", 
  "28 days"
)

ratio.late.candidate.markers.long <- ratio.late.candidate.markers.long %>% 
  mutate(`wound age` = factor(`wound age`, levels = wound_age_order)) %>% 
  arrange(`wound age`) %>% 
  mutate(sample = factor(sample, levels = unique(sample)))


plots <- lapply(unique(ratio.late.candidate.markers.long$gene), function(x){
  ratio.late.candidate.markers.long %>% 
    filter(gene == x) %>%
    ggplot(aes(x = sample, y = ratio)) +
    geom_point(size = 5, color = "steelblue") +
    geom_text(
      aes(label = paste(localisation, "\n", `Wound age`, age, Gender, sep = ", ")),
      vjust = -0.5, size = 2.5
    ) +
    labs(
      title = paste("Expression for", x),
      x = "",
      y = "Expression Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    #expand_limits(y = max(ratio) * 1.1) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 14),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 20) 
    ) + 
    coord_cartesian(clip = "off")
})

for (i in seq_along(plots)) {
  ggsave(file.path(out_data, paste0("04_plot_", unique(ratio.late.candidate.markers.long$gene)[i], ".late.pdf")), plots[[i]], width = 8, height = 6)
}


## Final heatmap with all candidates. Plot ratios 
genes <- c(row.names(ratio.early.markers),
           row.names(ratio.int.markers),
           row.names(ratio.late.markers)
)

sample_order <- c(levels(ratio.early.candidate.markers.long$sample),levels(ratio.int.candidate.markers.long$sample), levels(ratio.late.candidate.markers.long$sample))

ratio.candiates <- ratio[genes, sample_order]

library(ComplexHeatmap)
library(circlize) # for color functions

mat <- ratio.candiates
mat[mat < 1] <- -1 / mat[mat < 1]

mat[mat > 20] <- 20
mat[mat < -20] <- -20


stage_anno <- HeatmapAnnotation(
  Stage = factor(
    c(rep("Early", length(levels(ratio.early.candidate.markers.long$sample))),
      rep("Intermediary", length(levels(ratio.int.candidate.markers.long$sample))),
      rep("Late", length(levels(ratio.late.candidate.markers.long$sample)))),
    levels = c("Early", "Intermediary", "Late")
  ),
  col = list(Stage = c("Early" = "#4daf4a", 
                       "Intermediary" = "#377eb8", 
                       "Late" = "#e41a1c"))
)

col_fun <- colorRamp2(
  c(min(mat), 0, max(mat)),
  c("royalblue", "white", "tomato1")
)


# Replace with your actual gene sets
early_genes        <- row.names(ratio.early.markers) # replace with your early markers
intermediary_genes <- row.names(ratio.int.markers)  # replace with your intermediary markers
late_genes         <- row.names(ratio.late.markers)  # replace with your late markers


# Assign group to each gene in the same order as in mat
gene_stage <- ifelse(rownames(mat) %in% early_genes,        "Early",
                     ifelse(rownames(mat) %in% intermediary_genes, "Intermediary",
                            ifelse(rownames(mat) %in% late_genes,         "Late", NA)))


row_anno <- rowAnnotation(
  Marker = factor(gene_stage, levels = c("Early", "Intermediary", "Late")),
  col = list(Marker = c(
    "Early"        = "#4daf4a", 
    "Intermediary" = "#377eb8", 
    "Late"         = "#e41a1c"
  ))
)

pdf(file.path(out_data, "05_marker_heatmap.pdf"), w = 8, h =10)
Heatmap(
  mat,
  name = "Ratio",
  top_annotation = stage_anno,
  left_annotation = row_anno,  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = sample_order,
  col = col_fun
)
dev.off()
