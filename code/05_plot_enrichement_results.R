#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Plot Enrichment results from ORSUM and other ontologies
# @software version: R=4.4.2

# Library
sh <- suppressPackageStartupMessages

sh(library(tidyverse))
sh(library(DOSE))

out.dir.results <- "results_new/05_enrichement_plots/"

if (!dir.exists(out.dir.results)){
  dir.create(out.dir.results)
}

load_and_process_orsum <- function(filepath, comparison_name) {
  df <- read.csv(filepath, sep = "\t")
  colnames(df)[6] <- "Comparison Rank"
  df %>% mutate(comparison = comparison_name)
  #df %>% slice_head(n = 15) 
}

# Define a helper function to simplify loading with direction label
load_orsum <- function(path, contrast, direction) {
  load_and_process_orsum(path, contrast) %>% mutate(direction = direction)
}

# Analysis of the Up + down ----
# (where Up and down are separated) 
## Load ORSUM enrichment results ----
orsum_terms <- do.call("rbind", list(
  load_orsum("results_new/03_enrichement/Control - Early/ORSUM_bp_up/filteredResult-Summary.tsv", 
             "Control vs Early", "Up"),
  
  load_orsum("results_new/03_enrichement/Control - Intermediary/ORSUM_bp_up/filteredResult-Summary.tsv", 
             "Control vs Intermediary", "Up"),
  
  load_orsum("results_new/03_enrichement/Control - Intermediary/ORSUM_bp_down/filteredResult-Summary.tsv", 
             "Control vs Intermediary", "Down"),
  
  load_orsum("results_new/03_enrichement/Control - Late/ORSUM_bp_up/filteredResult-Summary.tsv", 
             "Control vs Late", "Up"),
  
  load_orsum("results_new/03_enrichement/Control - Late/ORSUM_bp_down/filteredResult-Summary.tsv", 
             "Control vs Late", "Down"),
  
  load_orsum("results_new/03_enrichement/Early vs Intermediary/ORSUM_bp_down/filteredResult-Summary.tsv", 
             "Early vs Intermediary", "Down"),
  
  load_orsum("results_new/03_enrichement/Intermediary vs Late/ORSUM_bp_up/filteredResult-Summary.tsv", 
             "Intermediary vs Late", "Up"),
  
  load_orsum("results_new/03_enrichement/Intermediary vs Late/ORSUM_bp_down/filteredResult-Summary.tsv", 
             "Intermediary vs Late", "Down")
))


# Filter the most important terms per analysis, keep only those with Ranking <= 15
orsum_terms <- orsum_terms %>% filter(Representing.term.rank < 15)


#The issue with this approach is the fact that some terms that are the same across enrichement have different rank and might disapear 
#Instead lets jointly summarise 
orsum_terms <- read.csv("results_new/03_enrichement/Summary_All_Contrasts/filteredResult-Summary.tsv", sep = "\t") %>% 
  filter(Representing.term.rank < 11) %>% 
  rename(repr_rank = `Representing.term.rank`)

orsum_terms <- orsum_terms %>%
  pivot_longer(
    cols = ends_with(".term.rank"),  # Select all rank columns
    names_to = "contrast",
    values_to = "rank"
  )
orsum_terms <- orsum_terms %>% filter(rank != "None")


# Add direction
orsum_terms$direction <- ifelse(grepl("^Up", orsum_terms$contrast), "Up", "Down")

#Replace name in the contrast
contrast_vector <- c("Early vs Intermediary", "Control vs Intermediary", "Intermediary vs Late", "Control vs Late", 
                     "Control vs Early", "Control vs Intermediary", "Intermediary vs Late", "Control vs Late")
names(contrast_vector) <- names(table(orsum_terms$contrast))


orsum_terms$contrast <- unname(contrast_vector[orsum_terms$contrast])

# Order
orsum_terms$contrast <- factor(orsum_terms$contrast, levels = c("Control vs Early", "Control vs Intermediary",
                                                                "Early vs Intermediary", "Control vs Late", "Intermediary vs Late"))
orsum_terms$direction <- factor(orsum_terms$direction, levels = c("Up", "Down"))


orsum_terms <- orsum_terms %>% 
  arrange(contrast, direction)

orsum_terms$`Representing.term.name` <- factor(orsum_terms$`Representing.term.name`,
                                               levels = unique(orsum_terms$`Representing.term.name`)[length( unique(orsum_terms$`Representing.term.name`)):1])




## Plot (GO) Orsum data ----
pdf(paste0(out.dir.results, "GeneOntology_BP_Orsum.pdf"), w = 15, h = 12)
plot.go.bp.dir <- ggplot(orsum_terms, aes(x = contrast, y = `Representing.term.name`, colour = direction)) +
  geom_point(size = 4) +
  scale_colour_manual(name = "Gene Set", values =  c("tomato3", "royalblue")) +
  theme_dose(10) +
  ylab("") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        plot.margin = margin(l = 100), 
        legend.text = element_text(size = 16),
        legend.title = element_blank()) + 
  xlab("")
plot(plot.go.bp.dir)
dev.off()



## Load KEGG/DO/DNG data----
load_enrichement_data <- function(filepath, comparison_name){
  df <- read.csv(filepath, sep = ",")
  df %>% mutate(comparison = comparison_name) %>% slice_head(n = 10)
}

### KEGG ----
kegg_terms_early.up <- load_enrichement_data(
  "results_new/03_enrichement/Control - Early/kegg.up.results.csv",
  "Control vs Early"
) %>% mutate(direction = "Up")

kegg_terms_intemediary.up <- load_enrichement_data(
  "results_new/03_enrichement/Control - Intermediary/kegg.up.results.csv",
  "Control vs Intermediary"
) %>% mutate(direction = "Up")

kegg_terms_intemediary.down <- load_enrichement_data(
  "results_new/03_enrichement/Control - Intermediary/kegg.down.results.csv",
  "Control vs Intermediary"
) %>% mutate(direction = "Down")

kegg_terms_late.up <- load_enrichement_data(
  "results_new/03_enrichement/Control - Late/kegg.up.results.csv",
  "Control vs Late"
) %>% mutate(direction = "Up")

kegg_terms_early_intemediary.up <- load_enrichement_data(
  "results_new/03_enrichement/Early vs Intermediary/kegg.up.results.csv",
  "Early vs Intermediary"
) %>% mutate(direction = "Up")

kegg_terms_early_intemediary.down <- load_enrichement_data(
  "results_new/03_enrichement/Early vs Intermediary/kegg.down.results.csv",
  "Early vs Intermediary"
) %>% mutate(direction = "Down")

kegg_terms_intemediary_late.up <- load_enrichement_data(
  "results_new/03_enrichement/Intermediary vs Late/kegg.up.results.csv",
  "Intermediary vs Late"
) %>% mutate(direction = "Up")

kegg_terms_intemediary_late.down <- load_enrichement_data(
  "results_new/03_enrichement/Intermediary vs Late/kegg.down.results.csv",
  "Intermediary vs Late"
) %>% mutate(direction = "Down")

KEGG <- do.call("rbind", 
                list(kegg_terms_early.up,
                     kegg_terms_intemediary.up,
                     kegg_terms_intemediary.down,
                     kegg_terms_late.up,
                     kegg_terms_early_intemediary.up,
                     kegg_terms_early_intemediary.down,
                     kegg_terms_intemediary_late.up,
                     kegg_terms_intemediary_late.down
                )
)


KEGG <- KEGG %>% 
  mutate(comparison = factor(comparison, levels = c("Control vs Early", "Control vs Intermediary",
                                                     "Early vs Intermediary", "Control vs Late", "Intermediary vs Late")),
         direction = factor(direction, levels = c("Up", "Down"))) %>% 
  arrange(comparison,direction) %>% 
  mutate(Description = factor(Description, levels = unique(Description)[length(unique(Description)):1]))

plot.KEGG.dir <- ggplot(KEGG, aes(x = comparison, y = Description, colour = direction)) +
  geom_point(size = 4) +
  scale_colour_manual(name = "Gene Set", values =  c("tomato3", "royalblue")) +
  theme_dose(10) +
  ylab("") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        plot.margin = margin(l = 100), 
        legend.text = element_text(size = 16),
        legend.title = element_blank()) + 
  xlab("")

pdf(paste0(out.dir.results, "KEGG.direction.pdf"), w = 15, h = 12)
plot(plot.KEGG.dir)
dev.off()
